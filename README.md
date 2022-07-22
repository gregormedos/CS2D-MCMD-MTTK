# CS2D-MCMD-MTTK
Monte Carlo (MC) and molecular dynamics (MD) programs for simulating a Core-softened (CS) anomalous liquid with PBC in 2 dimensions (2D) using the Martyna-Tuckerman-Tobias-Klein (MTTK) thermostat and barostat.

Compile Fortran90 program with gfortran (other compilers have not been tested!):

$ gfortran -O3 -o cs2d_nvt_cutoff cs2d_nvt_cutoff.f90

or

$ gfortran -O3 -o cs2d_npt_cutoff cs2d_npt_cutoff.f90

Create folder with the temperature/pressure value for the name:

$ mkdir 1.000

Run Python script to start a batch of jobs for a series of densities/temperatures at the given temperature/pressure:

$ python3 cs2d_nvt_cutoff.py

or

$ python3 cs2d_npt_cutoff.py
