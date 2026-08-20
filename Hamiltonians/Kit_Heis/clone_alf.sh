#!/bin/sh
BRANCH="master"
HAMILTONIAN_NAME="Kit_Heis"

git clone --depth 1 https://github.com/ALF-QMC/ALF.git --branch $BRANCH || exit 1
ln -s "$PWD/Hamiltonian_${HAMILTONIAN_NAME}_smod.F90" ALF/Prog/Hamiltonians || exit 1
echo "$HAMILTONIAN_NAME" > "ALF/Prog/Hamiltonians.list.d/$HAMILTONIAN_NAME"
echo "Finished cloning ALF and setting up the Hamiltonian \"$HAMILTONIAN_NAME\"."

echo "Creating symbolic links and modifying the Makefile for the analysis that calculates the magetotropic susceptibility."
echo "The workflow that adapts the Makefile is fragile and may break if the Makefile is changed in future versions of ALF."
ln -s "$PWD/calc_k2_tau.F90" "ALF/Analysis/calc_k2_tau.F90" || exit 1
sed -i "s|BINS=|BINS= calc_k2_tau.out |" ALF/Analysis/Makefile || exit 1
echo "calc_k2_tau.o: ana_mod.o" >> ALF/Analysis/Makefile || exit 1
echo "Done."
