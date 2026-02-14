#!/bin/bash
# Run all hexabody registered-acquisition cells in parallel on descent.
# Usage: ssh descent "cd ~/julia_shared_dev/SMLMAnalysis && bash dev/run_hexabody_all.sh"
#
# Runs N_PARALLEL jobs at a time (default 2). Each cell takes ~15-45 min.
# GPU is shared across processes (detectfit serializes on CUDA).

set -euo pipefail
N_PARALLEL=${1:-2}
SCRIPT="dev/genmab_hexabody_cell.jl"
LOGDIR="dev/output/logs"
mkdir -p "$LOGDIR"

# Build job list: h5file outname
declare -a JOBS=()

DATA="/mnt/nas/cellpath/Genmab/Data"

# --- 20250603 ---
BASE="$DATA/20250603_A431_SaturatingIgG10min+C1q"

# 2F8 wild-type
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0603_wt_cell${i}")
done

# E345R
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-E345R-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0603_e345r_cell${i}")
done

# E430G
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-E430G-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0603_e430g_cell${i}")
done

# RGY
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-RGY-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0603_rgy_cell${i}")
done

# --- 20250611 ---
BASE="$DATA/20250611_A431_SaturatingIgG10min+C1q"

# 2F8 wild-type
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0611_wt_cell${i}")
done

# E345R (Cell_01 and Cell_02 have aborted small files - tail -1 picks the later/larger one)
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-E345R-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0611_e345r_cell${i}")
done

# E430G
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-E430G-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0611_e430g_cell${i}")
done

# RGY
for i in $(seq -w 1 11); do
    h5=$(ls "$BASE/A431_IgG1-2F8-RGY-AF647_5ugml_10min+C1q/Cell_${i}/Label_01/"*.h5 | tail -1)
    JOBS+=("$h5 0611_rgy_cell${i}")
done

echo "Total cells: ${#JOBS[@]}"
echo "Parallel jobs: $N_PARALLEL"
echo ""

# Run with limited parallelism
running=0
pids=()
for job in "${JOBS[@]}"; do
    read -r h5file outname <<< "$job"
    logfile="$LOGDIR/${outname}.log"

    # Skip if output already exists and has a final render
    if [ -d "dev/output/${outname}/07_render" ]; then
        echo "SKIP $outname (already complete)"
        continue
    fi

    # Skip aborted/tiny files (<1GB)
    fsize=$(stat -c%s "$h5file" 2>/dev/null || echo 0)
    if [ "$fsize" -lt 1073741824 ]; then
        echo "SKIP $outname (file too small: $(( fsize / 1048576 ))MB, likely aborted)"
        continue
    fi

    echo "START $outname -> $logfile"
    julia -t auto --project=dev "$SCRIPT" "$h5file" "$outname" > "$logfile" 2>&1 &
    pids+=($!)
    running=$((running + 1))

    # Wait if at capacity
    if [ "$running" -ge "$N_PARALLEL" ]; then
        wait -n  # wait for any one to finish
        running=$((running - 1))
    fi
done

# Wait for remaining
for pid in "${pids[@]}"; do
    wait "$pid" 2>/dev/null || true
done

echo ""
echo "All jobs complete. Check $LOGDIR/ for per-cell logs."

# Quick summary
echo ""
echo "=== RESULTS SUMMARY ==="
for job in "${JOBS[@]}"; do
    read -r _ outname <<< "$job"
    statsfile="dev/output/${outname}/05_driftcorrect/stats.md"
    if [ -f "$statsfile" ]; then
        drift=$(grep "Max intra-dataset drift" "$statsfile" | grep -oP '[\d.]+' || echo "?")
        shift=$(grep "Max inter-dataset shift" "$statsfile" | grep -oP '[\d.]+' || echo "?")
        printf "%-25s drift=%snm  shift=%snm\n" "$outname" "$drift" "$shift"
    else
        printf "%-25s FAILED\n" "$outname"
    fi
done
