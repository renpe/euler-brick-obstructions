"""
Dispatcher for the full torsion-intersection scan.

Generates all (m,n) coprime pairs with 2 <= m <= M_MAX, n < m, m - n odd,
splits them into NWORKERS chunks (round-robin), and starts NWORKERS Sage
subprocesses, each processing one chunk. Collects JSONL outputs and
aggregates them into a CSV of proven fibers.

Usage:
    python jh_torsion_full_dispatch.py [M_MAX] [NWORKERS]

Environment:
    SAGE_BIN  — path to the sage binary (default: /home/varon/miniforge3/envs/sage/bin/sage)
    OUT_DIR   — output directory (default: ./scan_output_M{M_MAX})
"""
import os
import sys
import csv
import json
import subprocess
import shutil
from math import gcd
from time import time

SAGE_BIN = os.environ.get(
    'SAGE_BIN',
    '/home/varon/miniforge3/envs/sage/bin/sage',
)

HERE = os.path.dirname(os.path.abspath(__file__))
WORKER = os.path.join(HERE, 'jh_torsion_full_worker.py')


def generate_pairs(m_max):
    pairs = []
    for m in range(2, m_max + 1):
        for n in range(1, m):
            if gcd(m, n) != 1:
                continue
            if (m - n) % 2 != 1:
                continue
            pairs.append((m, n))
    return pairs


def main():
    m_max = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    n_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 30

    out_dir = os.environ.get(
        'OUT_DIR',
        os.path.join(HERE, '..', 'data', f'scan_output_M{m_max}'),
    )
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    pairs = generate_pairs(m_max)
    print(f"[dispatch] M_MAX={m_max}, total coprime pairs: {len(pairs)}",
          flush=True)
    print(f"[dispatch] starting {n_workers} workers, output -> {out_dir}",
          flush=True)

    # Round-robin split
    chunks = [[] for _ in range(n_workers)]
    for i, p in enumerate(pairs):
        chunks[i % n_workers].append(p)

    # Write chunk files
    chunk_files = []
    out_files = []
    for w in range(n_workers):
        cf = os.path.join(out_dir, f'chunk_{w:02d}.csv')
        of = os.path.join(out_dir, f'out_{w:02d}.jsonl')
        with open(cf, 'w') as f:
            for (m, n) in chunks[w]:
                f.write(f'{m},{n}\n')
        chunk_files.append(cf)
        out_files.append(of)

    # Launch workers
    t0 = time()
    procs = []
    for w in range(n_workers):
        log_file = os.path.join(out_dir, f'worker_{w:02d}.log')
        log_fp = open(log_file, 'w')
        cmd = [SAGE_BIN, '-python', WORKER, chunk_files[w], out_files[w]]
        proc = subprocess.Popen(cmd, stdout=log_fp, stderr=log_fp)
        procs.append((proc, log_fp))
    print(f"[dispatch] all {n_workers} workers launched", flush=True)

    # Wait for all
    for w, (proc, log_fp) in enumerate(procs):
        rc = proc.wait()
        log_fp.close()
        if rc != 0:
            print(f"[dispatch] worker {w} exited rc={rc}", flush=True)

    elapsed = time() - t0
    print(f"[dispatch] all workers done after {elapsed:.0f}s "
          f"({elapsed/60:.1f}min)", flush=True)

    # Aggregate results
    all_results = []
    for of in out_files:
        if not os.path.exists(of):
            continue
        with open(of) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    all_results.append(json.loads(line))
                except json.JSONDecodeError:
                    continue

    # Write proven_fibers CSV (proven fibers only)
    csv_path = os.path.join(HERE, '..', 'data', f'proven_fibers_M{m_max}.csv')
    csv_path = os.path.abspath(csv_path)
    n_proven = 0
    by_quot = {}
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['m', 'n', 'quotient', 'method', 'tors_size',
                    'rk_method_quot', 'rk_PQ', 'rk_uV', 'rk_3'])
        for r in all_results:
            if not r.get('proven'):
                continue
            n_proven += 1
            quot = r.get('proven_via')
            method = r.get('proven_method', 'naive')
            by_quot[(quot, method)] = by_quot.get((quot, method), 0) + 1
            tors = r.get(f'tors_{quot}')
            rk_method = r.get(f'rk_method_{quot}')
            w.writerow([r['m'], r['n'], quot, method, tors, rk_method,
                        r.get('rk_PQ'), r.get('rk_uV'), r.get('rk_3')])

    full_path = os.path.join(HERE, '..', 'data', f'scan_full_M{m_max}.jsonl')
    full_path = os.path.abspath(full_path)
    with open(full_path, 'w') as f:
        for r in all_results:
            f.write(json.dumps(r) + '\n')

    print(f"\n========== SUMMARY ==========", flush=True)
    print(f"M_MAX:              {m_max}", flush=True)
    print(f"Total fibers:       {len(all_results)}", flush=True)
    print(f"Proven total:       {n_proven}", flush=True)
    print(f"  by quotient: {by_quot}", flush=True)
    print(f"Proven CSV:         {csv_path}", flush=True)
    print(f"Full JSONL:         {full_path}", flush=True)
    print(f"Elapsed:            {elapsed/60:.1f} min", flush=True)


if __name__ == '__main__':
    main()
