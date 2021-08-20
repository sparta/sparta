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
    exit 0
fi

echo "STATUS: Reblessing log files..."
if [ ! "$1" = "--rerun-failed" ]; then
    ################################################################################
    rm -rf CMake*
    ################################################################################
    cmake -C ../cmake/presets/mac_mpi.cmake \
        -DSPARTA_ENABLE_TESTING=ON \
        -DSPARTA_DSMC_TESTING_DRIVER_ARGS='-auto-rebless true' \
        ../cmake
    ################################################################################
    make -j4
    ################################################################################
    #ctest
else
    ctest --rerun-failed
fi

################################################################################
oldDateStr=24Aug20
dateStr=$(date "+%d%b%y")
for logFileOld in $(ls examples/**/log.$oldDateStr.*); do
    logFileNew=$(echo $logFileOld | sed "s/.${oldDateStr}//g")
    logFileGold=$(echo $logFileNew | sed "s/.mpi/.${dateStr}.mpi/g")
    #echo $logFileOld $logFileNew $logFileGold
    git rm ../$logFileOld
    mv $logFileNew ../$logFileGold
    git add ../$logFileGold
done
################################################################################
#git commit -m "examples: Reblessed log.archive.$dateStr"
################################################################################
echo "STATUS: Log files re-blessed. Please review 'git diff --staged'"
