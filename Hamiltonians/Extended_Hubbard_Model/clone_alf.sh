#!/bin/sh
BRANCH="master"
HAMILTONIAN_NAME="Extended_Hubbard"

git clone https://git.physik.uni-wuerzburg.de/ALF/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > ALF/Prog/Hamiltonians.list.d/Extended_Hubbard
