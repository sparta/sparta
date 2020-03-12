# Quick Start
```bash
cd /path/to/sparta
mkdir build
cd build
cmake ..
cmake -LH
cmake -LAH
make
```

# Quick start build triaging
``bash
cmake --log-level=VERBOSE [-C ../cmake/machine_files/FILE.cmake] ..
make VERBOSE=1
```
