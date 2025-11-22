#!/bin/bash
#
# Benchmark script for comparing old vs new svfits versions
#

# Paths to the two versions
NEW_SVFITS="/lustre_archive/apps/astrosoft/test_imaging/test_svfits/svfits"
OLD_SVFITS="/lustre_archive/apps/astrosoft/test_imaging/svfits/svfits"

# Parameter file (default or from command line)
PARFILE="${1:-svfits_par.txt}"

# Output files for each version
NEW_FITS="TEST_NEW.FITS"
OLD_FITS="TEST_OLD.FITS"

# Number of runs for averaging (default 3)
NRUNS="${2:-3}"

echo "=============================================="
echo "       SVFITS Benchmark Script"
echo "=============================================="
echo ""
echo "Parameter file: $PARFILE"
echo "Number of runs: $NRUNS"
echo ""

# Check if parameter file exists
if [ ! -f "$PARFILE" ]; then
    echo "ERROR: Parameter file '$PARFILE' not found!"
    exit 1
fi

# Check if both executables exist
if [ ! -x "$NEW_SVFITS" ]; then
    echo "ERROR: New svfits not found at $NEW_SVFITS"
    exit 1
fi

if [ ! -x "$OLD_SVFITS" ]; then
    echo "ERROR: Old svfits not found at $OLD_SVFITS"
    exit 1
fi

# Function to run benchmark
run_benchmark() {
    local SVFITS_BIN="$1"
    local OUTPUT_FITS="$2"
    local LABEL="$3"

    echo "----------------------------------------------"
    echo "Benchmarking: $LABEL"
    echo "Binary: $SVFITS_BIN"
    echo "----------------------------------------------"

    # Arrays to store results
    declare -a times
    declare -a memories
    declare -a cpu_usages

    for ((i=1; i<=NRUNS; i++)); do
        echo "  Run $i/$NRUNS..."

        # Remove old output
        rm -f TEST.FITS svfits.log 2>/dev/null

        # Run with /usr/bin/time to get memory and time stats
        # -v gives verbose output including max RSS
        OUTPUT=$(/usr/bin/time -v "$SVFITS_BIN" -u "$PARFILE" 2>&1)

        # Extract metrics from time output
        ELAPSED=$(echo "$OUTPUT" | grep "Elapsed (wall clock)" | awk '{print $NF}')
        MAX_RSS=$(echo "$OUTPUT" | grep "Maximum resident set size" | awk '{print $NF}')
        CPU_PERCENT=$(echo "$OUTPUT" | grep "Percent of CPU" | awk '{print $NF}' | tr -d '%')

        # Convert elapsed time to seconds (handles m:ss.ss and h:mm:ss formats)
        if [[ "$ELAPSED" == *":"* ]]; then
            # Has minutes or hours
            if [[ "$ELAPSED" == *":"*":"* ]]; then
                # h:mm:ss format
                SECS=$(echo "$ELAPSED" | awk -F: '{print ($1*3600)+($2*60)+$3}')
            else
                # m:ss format
                SECS=$(echo "$ELAPSED" | awk -F: '{print ($1*60)+$2}')
            fi
        else
            SECS="$ELAPSED"
        fi

        times+=("$SECS")
        memories+=("$MAX_RSS")
        cpu_usages+=("$CPU_PERCENT")

        echo "    Time: ${ELAPSED} (${SECS}s), Memory: ${MAX_RSS} KB, CPU: ${CPU_PERCENT}%"

        # Save output file from last run
        if [ -f "TEST.FITS" ] && [ $i -eq $NRUNS ]; then
            mv TEST.FITS "$OUTPUT_FITS"
        fi
    done

    # Calculate averages
    local sum_time=0
    local sum_mem=0
    local sum_cpu=0
    for t in "${times[@]}"; do sum_time=$(echo "$sum_time + $t" | bc -l); done
    for m in "${memories[@]}"; do sum_mem=$((sum_mem + m)); done
    for c in "${cpu_usages[@]}"; do sum_cpu=$(echo "$sum_cpu + $c" | bc -l); done

    AVG_TIME=$(echo "scale=3; $sum_time / $NRUNS" | bc -l)
    AVG_MEM=$((sum_mem / NRUNS))
    AVG_CPU=$(echo "scale=1; $sum_cpu / $NRUNS" | bc -l)

    # Store in global variables for comparison
    if [ "$LABEL" == "NEW (parallel)" ]; then
        NEW_AVG_TIME="$AVG_TIME"
        NEW_AVG_MEM="$AVG_MEM"
        NEW_AVG_CPU="$AVG_CPU"
    else
        OLD_AVG_TIME="$AVG_TIME"
        OLD_AVG_MEM="$AVG_MEM"
        OLD_AVG_CPU="$AVG_CPU"
    fi

    echo ""
    echo "  AVERAGE over $NRUNS runs:"
    echo "    Time:   ${AVG_TIME} seconds"
    echo "    Memory: ${AVG_MEM} KB ($(echo "scale=2; $AVG_MEM/1024" | bc) MB)"
    echo "    CPU:    ${AVG_CPU}%"
    echo ""
}

# Run benchmarks
run_benchmark "$NEW_SVFITS" "$NEW_FITS" "NEW (parallel)"
run_benchmark "$OLD_SVFITS" "$OLD_FITS" "OLD (serial)"

# Summary comparison
echo "=============================================="
echo "           SUMMARY COMPARISON"
echo "=============================================="
echo ""
printf "%-20s %15s %15s\n" "Metric" "NEW (parallel)" "OLD (serial)"
printf "%-20s %15s %15s\n" "------" "--------------" "-----------"
printf "%-20s %15.3f %15.3f\n" "Time (seconds)" "$NEW_AVG_TIME" "$OLD_AVG_TIME"
printf "%-20s %15d %15d\n" "Memory (KB)" "$NEW_AVG_MEM" "$OLD_AVG_MEM"
printf "%-20s %15.1f %15.1f\n" "CPU Usage (%)" "$NEW_AVG_CPU" "$OLD_AVG_CPU"
echo ""

# Calculate speedup and memory overhead
SPEEDUP=$(echo "scale=2; $OLD_AVG_TIME / $NEW_AVG_TIME" | bc -l)
MEM_RATIO=$(echo "scale=2; $NEW_AVG_MEM / $OLD_AVG_MEM" | bc -l)

echo "Speedup:          ${SPEEDUP}x faster"
echo "Memory ratio:     ${MEM_RATIO}x (new/old)"
echo ""

# Compare output files
echo "----------------------------------------------"
echo "Output File Comparison"
echo "----------------------------------------------"
if [ -f "$NEW_FITS" ] && [ -f "$OLD_FITS" ]; then
    NEW_SIZE=$(stat -c%s "$NEW_FITS" 2>/dev/null || stat -f%z "$NEW_FITS")
    OLD_SIZE=$(stat -c%s "$OLD_FITS" 2>/dev/null || stat -f%z "$OLD_FITS")

    echo "NEW FITS size: $NEW_SIZE bytes ($(echo "scale=2; $NEW_SIZE/1048576" | bc) MB)"
    echo "OLD FITS size: $OLD_SIZE bytes ($(echo "scale=2; $OLD_SIZE/1048576" | bc) MB)"

    if [ "$NEW_SIZE" -eq "$OLD_SIZE" ]; then
        echo "File sizes: MATCH"

        # Binary comparison
        if cmp -s "$NEW_FITS" "$OLD_FITS"; then
            echo "Binary comparison: IDENTICAL"
        else
            echo "Binary comparison: DIFFERENT (expected due to processing order)"
            # Try comparing with fitsdiff if available
            if command -v fitsdiff &> /dev/null; then
                echo ""
                echo "Running fitsdiff..."
                fitsdiff "$NEW_FITS" "$OLD_FITS"
            fi
        fi
    else
        echo "File sizes: DIFFERENT"
        echo "  Difference: $((NEW_SIZE - OLD_SIZE)) bytes"
    fi
else
    echo "WARNING: Could not compare output files"
    [ ! -f "$NEW_FITS" ] && echo "  Missing: $NEW_FITS"
    [ ! -f "$OLD_FITS" ] && echo "  Missing: $OLD_FITS"
fi

echo ""
echo "=============================================="
echo "Benchmark complete!"
echo "=============================================="
