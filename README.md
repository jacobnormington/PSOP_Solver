# PSOP_Solver
Parallel B&amp;B solver for the sequential ordering problem (SOP)

## Usage
- `make`.
- `./sop_solver <Instance Location> <Thread number> <Config_file>`
- Example: `./sop_solver soplib/R.500.1000.30.sop 4 soplib_config.txt`
The boost library is required before compilation. It can be installed via `sudo apt-get install libboost-all-dev`

The dynamic hungarian library code is referenced from https://github.com/rod409/SOP.
