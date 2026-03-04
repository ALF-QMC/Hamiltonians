#!/bin/sh
BRANCH="master"
HAMILTONIAN_NAME="Hubbard_Can"

git clone https://github.com/ALF-QMC/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > ALF/Prog/Hamiltonians.list.d/Hubbard_Can

