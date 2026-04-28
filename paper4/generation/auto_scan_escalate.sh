#!/bin/bash
# Phasen-Eskalation mit disjunkten m-Bereichen + Auto-Rerun von ★-Treffern.
#
# Phase 1: m ∈ [2, 200]
# Phase 2: m ∈ [201, 500]
# Phase 3: m ∈ [501, 1000]
#
# Nach jeder Phase: merge CSVs, extrahiere Rang≥2-Fasern, mw_rerun.py darauf.

set -u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate sage

# Run from this script's own directory (scripts/paper4/generation/)
cd "$(dirname "$0")"

LOG=/tmp/scan.log
NUM_WORKERS=4
PYTHON="${PYTHON:-../../../.venv-linux/bin/python}"

wait_for_phase_workers() {
    local M_MAX=$1
    while pgrep -f "search_new_fibres.sage $M_MAX " > /dev/null 2>&1; do
        sleep 30
    done
}

run_phase_workers() {
    local M_MIN=$1
    local M_MAX=$2
    echo "" >> "$LOG"
    echo "=== Phase m∈[$M_MIN, $M_MAX] gestartet $(date -Iseconds) ===" >> "$LOG"
    local pids=()
    for i in $(seq 0 $((NUM_WORKERS-1))); do
        sage search_new_fibres.sage $M_MAX $i $NUM_WORKERS $M_MIN \
            > /tmp/scan_w${i}_max${M_MAX}.log 2>&1 &
        pids+=($!)
        echo "  Worker $i (m%%${NUM_WORKERS}==$i, m∈[$M_MIN,$M_MAX]) PID=${pids[$i]}" >> "$LOG"
    done
    for pid in "${pids[@]}"; do
        wait $pid
    done
}

merge_and_rerun() {
    local M_MAX=$1
    local merged="rank_scan_max${M_MAX}.csv"
    echo "m,n,analytic_rank" > "$merged"
    for i in $(seq 0 $((NUM_WORKERS-1))); do
        local csv="rank_scan_w${i}of${NUM_WORKERS}_max${M_MAX}.csv"
        if [ -f "$csv" ]; then
            tail -n +2 "$csv" >> "$merged"
        fi
    done
    local hits=$(awk -F, 'NR>1 && $3>=2' "$merged" | wc -l)
    echo "" >> "$LOG"
    echo "=== Phase m≤$M_MAX merge fertig $(date -Iseconds), Rang≥2: $hits ===" >> "$LOG"
    awk -F, 'NR>1 && $3>=2 {print "  ★ ("$1","$2") rank="$3}' "$merged" >> "$LOG"

    if [ "$hits" -gt 0 ]; then
        local stars=$(awk -F, 'NR>1 && $3>=2 {print "["$1","$2"]"}' "$merged" | paste -sd ',')
        local stars_json="[$stars]"
        echo "  Auto-Rerun startet: $stars_json" >> "$LOG"
        # 4 Worker, 600s timeout, descent_limit=20, MAX_SCALAR=5 (mehr Punkte pro Faser)
        "$PYTHON" mw_rerun.py "$stars_json" 600 20 5 >> "$LOG" 2>&1
        echo "  Auto-Rerun fertig $(date -Iseconds)" >> "$LOG"
    fi
}

echo "" >> "$LOG"
echo "===== auto_scan_escalate_v2 gestartet $(date -Iseconds), Workers=$NUM_WORKERS =====" >> "$LOG"

# Phase 1 (m≤200): start fresh (Resume from CSV)
run_phase_workers 2 200
merge_and_rerun 200

# Phase 2: m∈[201, 500]
run_phase_workers 201 500
merge_and_rerun 500

# Phase 3: m∈[501, 1000]
run_phase_workers 501 1000
merge_and_rerun 1000

echo "" >> "$LOG"
echo "===== ALLE PHASEN FERTIG $(date -Iseconds) =====" >> "$LOG"
