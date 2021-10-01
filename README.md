`sesh` calculates corrections for time-of-flight experimental data.

------------------------
# Getting Started
------------------------

<u>Requirements</u>:
- Fortran compiler: `gfortran` (gcc 10.2)

<u>To compile</u>:
- Run `compile-sesh.sh`

<u>To run executable</u>:

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
fully defined in the report by F. Frohner described [here](https://www.osti.gov/biblio/4554018).

