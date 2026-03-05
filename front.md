# Repository of ALF Hamiltonians

Repository for Hamiltonians that are not part of the [ALF repository](https://github.com/ALF-QMC/ALF).

The ALF library comes with a number of supported model classes [^1] — this repository contains additional Hamiltonians, which may, if indicated, require a specific version of ALF.

[^1]: Currently ALF includes (i) SU(N) Hubbard models, (ii) O(2N) t-V models, (iii) Kondo models, (iv) long-range Coulomb models, and (v) generic Z2 lattice gauge theories coupled to Z2 matter and fermions. See its [documentation](https://alf.physik.uni-wuerzburg.de/page/documentation/) for descriptions of those models.


## List of Hamiltonians

```{tableofcontents}
```


## Running Hamiltonians

The  shell script will clone the ALF version that is required for the Hamiltonian and  will carry out the required changes, to the ALF code to account for the new Hamiltonian.

1. Clone this git project
2. Go to the directory corresponding to the desired model 
3. Execute the shell script `clone_alf.sh`[^2]
4. In the ALF  directory carry out the compilation as for a normal ALF implementation
5. The executable `$PWD/ALF/Prog/ALF.out` is enabled with the new Hamiltonian
6. The `Start` directory contains the initial files required to carry out a sample run

[^2]: Note that this creates a shallow clone of the `ALF` repository, meaning its full git history will not be fetched. If you want to access the whole history, run `git fetch --unshallow` inside the `ALF/` directory.

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
