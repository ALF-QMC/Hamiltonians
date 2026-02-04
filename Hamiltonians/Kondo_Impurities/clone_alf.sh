#!/bin/sh
BRANCH="Base_ALF_for_Kondo_Impurities"
HAMILTONIAN_NAME="Kondo_impurities"

git clone https://git.physik.uni-wuerzburg.de/ALF/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > ALF/Prog/Hamiltonians.list.d/Kondo_impurities

