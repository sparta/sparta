#!/bin/bash
if [ ! $(basename $(pwd)) = "build" ]; then
    echo "ERROR: This script must be run from the build directory"
    exit 1
fi

if [[ "$1" = "-h" || "$1" = "--help" ]]; then
    echo "Usage:"
    echo "  cd /path/to/sparta/build"
    echo "  /path/to/sparta/tools/testing/rebless.sh [--help] [--rerun-failed]"
    echo "Options:"
    echo "  --help: Print this help menu and exit"
    echo "  --rerun-failed: Use the existing build directory and rerun failed tests"
    echo "Environment:"
    echo "  SPARTA_REBLESS_PRESET: cmake preset to configure with"
    echo "                         (default: mpi)"
    exit 0
fi

preset=${SPARTA_REBLESS_PRESET:-mpi}

echo "STATUS: Reblessing log files..."
if [ ! "$1" = "--rerun-failed" ]; then
    ################################################################################
    rm -rf CMake*
    ################################################################################
    cmake -C ../cmake/presets/${preset}.cmake \
        -DSPARTA_ENABLE_TESTING=ON \
        -DSPARTA_DSMC_TESTING_DRIVER_ARGS='-auto-rebless true' \
        ../cmake
    ################################################################################
    make -j4
    ################################################################################
    ctest
else
    ctest --rerun-failed
fi

################################################################################
# Gold logs are named log.<date>.mpi_<ranks>.<problem>, and the date is the
# day that log was blessed, so it differs from one file to the next.  Take
# the date from each filename rather than assuming a single one for all of
# them.  A test run writes log.mpi_<ranks>.<problem> alongside its gold log;
# that is the file being promoted.
dateStr=$(date "+%d%b%y")
count=0
for logFileOld in $(ls examples/*/log.*.mpi_*.* 2>/dev/null); do
    dir=$(dirname $logFileOld)
    base=$(basename $logFileOld)
    # log . <date> . mpi_<ranks> . <problem>
    rest=${base#log.*.}
    logFileNew=$dir/log.$rest
    logFileGold=$dir/log.$dateStr.$rest
    if [ ! -f "$logFileNew" ]; then
        continue                      # this test did not run, keep its gold log
    fi
    git rm -q ../$logFileOld
    mv $logFileNew ../$logFileGold
    git add ../$logFileGold
    count=$((count+1))
done
################################################################################
#git commit -m "examples: Reblessed log.archive.$dateStr"
################################################################################
if [ $count -eq 0 ]; then
    echo "STATUS: No log files were re-blessed -- no test produced a new log."
    echo "        Check that the tests actually ran."
    exit 1
fi
echo "STATUS: $count log files re-blessed. Please review 'git diff --staged'"
