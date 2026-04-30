"""
Dispatcher for the mwrank pass.

Reads the existing scan_full_M{M_MAX}.jsonl, identifies fibers that
were not proven and have at least one quotient with rk=[0, k>0]
(ambiguous lower bound), and runs mwrank on each ambiguous quotient
in parallel via NWORKERS Sage subprocesses.

Usage:
    python jh_mwrank_pass_dispatch.py [M_MAX] [NWORKERS]

Environment:
    SAGE_BIN  — path to sage (default: /home/varon/miniforge3/envs/sage/bin/sage)
"""
import os
import sys
import json
import csv
import subprocess
from time import time

SAGE_BIN = os.environ.get(
    'SAGE_BIN', '/home/varon/miniforge3/envs/sage/bin/sage')

HERE = os.path.dirname(os.path.abspath(__file__))
WORKER = os.path.join(HERE, 'jh_mwrank_pass_worker.py')


def find_ambiguous(scan_jsonl_path):
    """Return list of (m, n, quotient) triples that need mwrank retesting."""
    triples = []
    with open(scan_jsonl_path) as f:
        for line in f:
            r = json.loads(line)
            if r.get('proven'):
                continue
            for q in ['E_3', 'E_uV', 'E_PQ']:
                rk = r.get(f'rk_{q}')
                if rk is None:
                    continue
                lo, hi = rk[0], rk[1]
                if lo == 0 and hi is not None and hi > 0:
                    triples.append((r['m'], r['n'], q))
    return triples


def main():
    m_max = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 30

    scan_jsonl = os.path.join(
        HERE, '..', 'data', f'scan_full_M{m_max}.jsonl')
    scan_jsonl = os.path.abspath(scan_jsonl)

    out_dir = os.path.join(
        HERE, '..', 'data', f'mwrank_pass_M{m_max}')
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    triples = find_ambiguous(scan_jsonl)
    print(f"[mwrank-pass] M_MAX={m_max}", flush=True)
    print(f"[mwrank-pass] ambiguous (m,n,quotient) triples: {len(triples)}",
          flush=True)
    if not triples:
        print("[mwrank-pass] nothing to do", flush=True)
        return

    chunks = [[] for _ in range(n_workers)]
    for i, t in enumerate(triples):
        chunks[i % n_workers].append(t)

    chunk_files = []
    out_files = []
    for w in range(n_workers):
        cf = os.path.join(out_dir, f'chunk_{w:02d}.csv')
        of = os.path.join(out_dir, f'out_{w:02d}.jsonl')
        with open(cf, 'w') as f:
            for (m, n, q) in chunks[w]:
                f.write(f'{m},{n},{q}\n')
        chunk_files.append(cf)
        out_files.append(of)

    print(f"[mwrank-pass] launching {n_workers} workers", flush=True)
    t0 = time()
    procs = []
    for w in range(n_workers):
        log_file = os.path.join(out_dir, f'worker_{w:02d}.log')
        log_fp = open(log_file, 'w')
        cmd = [SAGE_BIN, '-python', WORKER, chunk_files[w], out_files[w]]
        proc = subprocess.Popen(cmd, stdout=log_fp, stderr=log_fp)
        procs.append((proc, log_fp))

    for w, (proc, log_fp) in enumerate(procs):
        rc = proc.wait()
        log_fp.close()
        if rc != 0:
            print(f"[mwrank-pass] worker {w} exited rc={rc}", flush=True)

    elapsed = time() - t0
    print(f"[mwrank-pass] all workers done after {elapsed:.0f}s "
          f"({elapsed/60:.1f}min)", flush=True)

    # Aggregate
    mwrank_results = {}
    for of in out_files:
        if not os.path.exists(of):
            continue
        with open(of) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                key = (r['m'], r['n'], r['quotient'])
                mwrank_results[key] = r

    # Categorize
    n_certain_zero = 0
    n_certain_pos = 0
    n_uncertain = 0
    n_error = 0
    for r in mwrank_results.values():
        if 'error' in r:
            n_error += 1
        elif r.get('mwrank_certain'):
            if r['mwrank_rank'] == 0:
                n_certain_zero += 1
            else:
                n_certain_pos += 1
        else:
            n_uncertain += 1

    # Now merge: re-load original scan, apply mwrank results, and re-evaluate proven status.
    with open(scan_jsonl) as f:
        all_results = [json.loads(line) for line in f if line.strip()]

    proved_now = 0
    for r in all_results:
        if r.get('proven'):
            continue
        for q in ['E_3', 'E_uV', 'E_PQ']:
            key = (r['m'], r['n'], q)
            mw = mwrank_results.get(key)
            if not mw or 'error' in mw or not mw.get('mwrank_certain'):
                continue
            if mw['mwrank_rank'] == 0 and mw['mwrank_rank_bound'] == 0:
                # mwrank certified rank=0; update rank info
                r[f'rk_{q}'] = [0, 0]
                r[f'rk_method_{q}'] = 'mwrank'
                # If torsion already known, use it; else we cannot resolve
                tors = r.get(f'tors_{q}')
                if tors is not None:
                    target = {'E_PQ': 6, 'E_uV': 6, 'E_3': 4}[q]
                    if tors == target:
                        r['proven'] = True
                        r['proven_via'] = q
                        r['proven_method'] = 'mwrank_naive'
                        proved_now += 1
                        break

    # Write merged results
    merged_jsonl = os.path.join(
        HERE, '..', 'data', f'scan_full_M{m_max}_with_mwrank.jsonl')
    merged_jsonl = os.path.abspath(merged_jsonl)
    with open(merged_jsonl, 'w') as f:
        for r in all_results:
            f.write(json.dumps(r) + '\n')

    n_total_proven = sum(1 for r in all_results if r.get('proven'))

    print(f"\n========== MWRANK PASS SUMMARY ==========", flush=True)
    print(f"Triples processed:          {len(triples)}", flush=True)
    print(f"  certain rank=0:           {n_certain_zero}", flush=True)
    print(f"  certain rank>=1:          {n_certain_pos}", flush=True)
    print(f"  uncertain (still ambig):  {n_uncertain}", flush=True)
    print(f"  errors/timeouts:          {n_error}", flush=True)
    print(f"\nNew fibers proven (after mwrank): {proved_now}", flush=True)
    print(f"Total proven now: {n_total_proven} / {len(all_results)} = "
          f"{100.0*n_total_proven/len(all_results):.1f}%", flush=True)
    print(f"\nMerged JSONL: {merged_jsonl}", flush=True)
    print(f"Wallclock:    {elapsed/60:.1f} min", flush=True)


if __name__ == '__main__':
    main()
