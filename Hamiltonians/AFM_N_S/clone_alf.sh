#!/bin/sh
BRANCH="300-hdf5-fix-libdir"
HAMILTONIAN_NAME="AFM_N_S"

git clone https://git.physik.uni-wuerzburg.de/ALF/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > ALF/Prog/Hamiltonians.list.d/Nematic_Dirac

