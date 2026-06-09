#!/bin/sh
BRANCH="Make_lattice-bigger-nnlist-option"
HAMILTONIAN_NAME="Kit_Heis"

git clone --depth 1 https://github.com/ALF-QMC/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > "ALF/Prog/Hamiltonians.list.d/$HAMILTONIAN_NAME"
