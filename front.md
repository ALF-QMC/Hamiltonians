# Repository of ALF Hamiltonians

Repository for Hamiltonians that are not part of the [ALF repository](https://git.physik.uni-wuerzburg.de/ALF/ALF).

The ALF library comes with a number of supported model classes [1] — this repository contains additional Hamiltonians, which may, if indicated, require an especific version of ALF.

[1] Currently ALF includes (i) SU(N) Hubbard models, (ii) O(2N) t-V models, (iii) Kondo models, (iv) long-range Coulomb models, and (v) generic Z2 lattice gauge theories coupled to Z2 matter and fermions. See its [documentation](https://gitpages.physik.uni-wuerzburg.de/ALF/ALF_Webpage/page/documentation/) for descriptions of those models.


## List of Hamiltonians

```{tableofcontents}
```


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
