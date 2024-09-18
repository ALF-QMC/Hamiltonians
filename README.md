# Hamiltonians

Repository for Hamiltonians that are not part of the [ALF repository](https://git.physik.uni-wuerzburg.de/ALF/ALF), with automatically generated [website](http://gitpages.physik.uni-wuerzburg.de/alf/hamiltonians).

The website is generated using [Jupyter Book](https://jupyterbook.org/), with one documentation file per Hamiltonian. The documentation can be written in different flavors of Markdown (`README.md`), as Jupyter notebook (`README.ipynb`), or reStructuredText (`README.rst`).

## Running Hamiltonians

The  shell script will clone the ALF version that is required for the Hamiltonian and  will carry out the required changes, to the ALF code to account for the new Hamiltonian.

1. Clone this git project
2. Go to the directory corresponding to the desired model 
3. Execute the shell script  `clone_alf.sh`
4. In the ALF  directory carry out the compilation as for a normal ALF implementation
5. The executable `$PWD/ALF/Prog/ALF.out` is enabled with the new Hamiltonian
6. The `Start` directory contains the initial files required to carry out a sample run

## Directory structure

One directory per Hamiltonian:

```
Hamiltonians/
├── abc/
│   ├── README[.md|.ipynb|.rst]
│   └── Hamiltonian_abc_smod.F90
├── def/
|   ├── README[.md|.ipynb|.rst]
|   └── Hamiltonian_def_smod.F90
...
```

