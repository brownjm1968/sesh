`sesh` calculates corrections for time-of-flight experimental data.

The code `sesh` presented here is a modernized version of SESH re-written in Modern Fortran.
Original code is by F.H. Froehner published in the report: "SESH - A Fortran IV Code For Calculating 
the Self-shielding and Multiple Scattering Effects For Neutron Cross Section Data Interpretation in 
the Unresolved Resonance Region", Report GA-8380 (1968). This code has been adapted over the years 
from the input of several nuclear data researchers and their contributions are appreciated.

`sesh` solves the problem of resonance self-shielding, the experimental effect seen in time-of-flight
experiments where the measured time/energy-averaged transmission data cannot be compared directly to 
theoretical cross section due to the inequality: \<f(x)\> != f(\<x\>) for non-linear equations.

------------------------
# Getting Started
------------------------

**Requirements**:
- Fortran compiler: `gfortran` (gcc 10.2)

**To compile**:
- Run `compile-sesh.sh`

**To run executable**:

Open bash-like terminal and run:

```sh
sesh
```

which will prompt you for input and output file names. If you have already defined these names
in a dedicated file like `auto-in`, you can simply run:

```sh
sesh < auto-in
```

# Input Files

An input file called `ta181.inp` is given as an example for SESH input. The SESH input is more
fully defined in the report by F. Froehner described [here](https://www.osti.gov/biblio/4554018).

