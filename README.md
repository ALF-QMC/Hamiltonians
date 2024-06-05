# Hamiltonians

Repository for Hamiltonians that are not part of the [ALF repository](https://git.physik.uni-wuerzburg.de/ALF/ALF), with automatically generated [website](http://gitpages.physik.uni-wuerzburg.de/alf/hamiltonians).

The website is generated using [Jupyert Book](https://jupyterbook.org/), with one documentation file per Hamiltonian. The documentation can be written in different flavors of Markdown (`README.md`), as Jupyter notebook (`README.ipynb`), or reStructuredText (`README.rst`).

## Running Hamiltonians

1. Clone this git project
2. Goto the directory corresponding to the desired model
3. Execute the shell script  clone_alf.sh

The  shell script will clone the ALF version that is  required  for  the Hamiltonain and  will carry  out  the  required  changes,  to  the ALF code to account  for the new  Hamiltonian.

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

