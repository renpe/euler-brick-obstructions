#!/bin/bash
# auto_refill.sh — keep the generation pipeline busy when the user is away.
#
# Polls CPU usage every $POLL seconds. When total CPU of mw_rerun /
# search_new_fibres / mw_fibre_worker drops below $THRESHOLD, launches
# the next wave from a queue of progressively higher MAX_SCALAR runs.
#
# Each launched wave is a separate mw_rerun.py invocation with its own
# log file in /tmp/refill_*.log. The watchdog itself logs to
# /tmp/auto_refill.log.
#
# Usage:
#   nohup scripts/paper4/generation/auto_refill.sh > /tmp/auto_refill.log 2>&1 &

set -u
cd "$(dirname "$0")/../../.."  # project root
PROJECT=$(pwd)
PYTHON="$PROJECT/.venv-linux/bin/python"
RERUN="$PROJECT/scripts/paper4/generation/mw_rerun.py"
LOG=/tmp/auto_refill.log

THRESHOLD=1500   # if total CPU < THRESHOLD%, start next wave
POLL=300         # poll every 5 minutes
MAX_PARALLEL=4   # max concurrent mw_rerun jobs (besides watchdog itself)

# ---------- Snapshot fiber lists once ----------
TIER34_FILE=/tmp/tier34_fibers.json
TIER2_FILE=/tmp/tier2_fibers.json
ALLMW_FILE=/tmp/allmw_fibers.json

# Refresh fiber lists from DB (in case Phase-1 added new ones)
"$PYTHON" - <<'PYEOF' || exit 1
import psycopg, json
c = psycopg.connect("host=192.168.178.63 dbname=euler user=euler password=euler").cursor()
c.execute("""
    SELECT mh.m, mh.n, count(*) AS cnt FROM pub.master_hits mh
    JOIN pub.hit_groups hg ON hg.hit_id = mh.id
    JOIN pub.generator_groups g ON g.id = hg.group_id
    WHERE g.short_name LIKE 'MW-%'
    GROUP BY mh.m, mh.n
""")
rows = [(int(m), int(n), int(c)) for m, n, c in c.fetchall()]
allmw = [[m, n] for m, n, _ in rows]
tier2 = [[m, n] for m, n, c in rows if 200 <= c < 50000]
tier34 = [[m, n] for m, n, c in rows if c < 200]
import os
for path, data in [
    ("/tmp/allmw_fibers.json", allmw),
    ("/tmp/tier2_fibers.json", tier2),
    ("/tmp/tier34_fibers.json", tier34),
]:
    with open(path, "w") as f: json.dump(data, f)
print(f"Snapshot: {len(allmw)} all, {len(tier2)} tier2, {len(tier34)} tier3+4")
PYEOF

# ---------- Wave queue: (label, fiber_file, timeout, descent, max_scalar) ----------
QUEUE=(
    "tier34_max13:$TIER34_FILE:600:25:13"
    "tier2_max9:$TIER2_FILE:900:25:9"
    "tier34_max15:$TIER34_FILE:600:25:15"
    "tier2_max11:$TIER2_FILE:900:25:11"
    "tier34_max17:$TIER34_FILE:900:25:17"
    "tier2_max13:$TIER2_FILE:1200:25:13"
    "tier34_max19:$TIER34_FILE:900:30:19"
    "tier2_max15:$TIER2_FILE:1200:30:15"
    "allmw_max5_re:$ALLMW_FILE:1800:30:5"   # combined backfill
    "tier34_max21:$TIER34_FILE:1200:30:21"
    "tier2_max17:$TIER2_FILE:1500:30:17"
)

cpu_load() {
    ps -eo pcpu,cmd 2>/dev/null \
      | grep -E "mw_rerun|search_new|mw_fibre_worker|mw_dispatcher" \
      | grep -v grep \
      | awk '{sum+=$1} END {print int(sum)}'
}

n_active_rerun() {
    # Count only the *main* mw_rerun.py invocations (ppid=1 because launched
    # via nohup and detached from this shell). Pool children share the same
    # cmdline but have ppid != 1.
    ps -eo pid,ppid,cmd | awk '$2==1 && /mw_rerun/ {n++} END {print n+0}'
}

STATE_FILE=/tmp/auto_refill.state
if [ -f "$STATE_FILE" ]; then
    idx=$(cat "$STATE_FILE")
    echo "$(date -Iseconds) Resume from idx=$idx (state file $STATE_FILE)" >> "$LOG"
else
    idx=0
fi
echo "$(date -Iseconds) auto_refill started, queue length=${#QUEUE[@]}, threshold=${THRESHOLD}%, poll=${POLL}s, idx=$idx" >> "$LOG"
echo "$(date -Iseconds) project=$PROJECT  python=$PYTHON" >> "$LOG"

while true; do
    cpu=$(cpu_load)
    n_proc=$(n_active_rerun)
    echo "$(date -Iseconds) CPU=${cpu}%  active_reruns=${n_proc}" >> "$LOG"

    if [ "$cpu" -lt "$THRESHOLD" ] && [ "$n_proc" -lt "$MAX_PARALLEL" ] && [ "$idx" -lt "${#QUEUE[@]}" ]; then
        IFS=":" read -r label fiber_file timeout descent max_scalar <<< "${QUEUE[$idx]}"
        if [ ! -s "$fiber_file" ]; then
            echo "$(date -Iseconds) SKIP $label: fiber file $fiber_file empty/missing" >> "$LOG"
            idx=$((idx + 1))
            continue
        fi
        json=$(cat "$fiber_file")
        out="/tmp/refill_${label}.log"
        nohup "$PYTHON" "$RERUN" "$json" "$timeout" "$descent" "$max_scalar" \
            > "$out" 2>&1 &
        pid=$!
        echo "$(date -Iseconds) Started $label (PID $pid, MAX_SCALAR=$max_scalar, descent=$descent), log=$out" >> "$LOG"
        idx=$((idx + 1))
        echo "$idx" > "$STATE_FILE"
        sleep 60   # let new job ramp up before next CPU check
    fi

    if [ "$idx" -ge "${#QUEUE[@]}" ] && [ "$n_proc" -eq 0 ] && [ "$cpu" -lt 200 ]; then
        echo "$(date -Iseconds) Queue exhausted and pipeline idle. Exiting watchdog." >> "$LOG"
        exit 0
    fi

    sleep "$POLL"
done
