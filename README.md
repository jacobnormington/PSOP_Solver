# PSOP_Solver
Sequential B&amp;B solver for the sequential ordering problem (SOP)

## Usage
- `make`.
- `./sop_solver <Instance Location> <Time limit[s]>`
- Example: `./sop_solver soplib/R.500.1000.30.sop 3600`
The boost library is required before compilation. It can be installed via `sudo apt-get install libboost-all-dev`