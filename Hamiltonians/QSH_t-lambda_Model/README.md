# QSH t - λ model


## Model definition

The model is defined and discussed in [this article](https://arxiv.org/abs/2410.05059) {cite}`Rein24` and is implemented as 

$$
\hat{H}  =  -t\sum_{\langle\boldsymbol{i},\boldsymbol{j}\rangle}\sum_{\sigma=\uparrow,\downarrow}\sum_{s=1}^{N_{\textrm{SU}(N)}}\big(\hat{c}_{\boldsymbol{i},\sigma,s}^\dagger\hat{c}_{\boldsymbol{j},\sigma,s}+\textrm{h.c.}\big) 
-\frac{\lambda}{N_{\textrm{SU}(N)}}\sum_{\boldsymbol{r}_{\textrm{Hexagon}}}\Big(\sum_{\sigma,\sigma'}\sum_{s=1}^{N_{\textrm{SU}(N)}}\sum_{\langle\langle\boldsymbol{i},\boldsymbol{j}\rangle\rangle\in \boldsymbol{r}_{\textrm{Hexagon}}}i\nu_{\boldsymbol{i}\boldsymbol{j}}\hat{c}_{\boldsymbol{i},\sigma,s}^\dagger\boldsymbol{\sigma}_{\sigma,\sigma'}\hat{c}_{\boldsymbol{j},\sigma',s}+\textrm{h.c.}\Big)^2
$$

Physically, the implementation supports the honeycomb lattice. Note that the spin is encoded as layer such that under `&VAR_lattice` in the parameter file `Lattice_type = "Bilayer_honeycomb"` has to be chosen.    The parameter file  for  this
specific model reads:

```fortran
&VAR_QSH              !! Variables for the specific model
ham_T      = 1.d0           ! Hopping parameter
ham_T2     = 1.d0
ham_chem   = 0.d0           ! Chemical potential
ham_lambda = 0.1d0
/
```
In the above Ham_T is  the nearest  neighbor  hopping and Ham_lambda is the coupling strength.  Finally Ham_chem  is  the chemical potential.
To use  this  Hamiltonian  you  have to specify:
```fortran
&VAR_ham_name
ham_name = "QSH"
/
```
in the parameters  file.

## Interaction
In order to reduce the Trotter error, the interaction term is implemented by using a Trotter decomposition corresponding to a Kekulé pattern by splitting up the lattice into three subgroups A,B,C. The specific implementation of the interaction is visually sketched in the `smod.F90` file of the Hamiltonian.

## Observables 
The code  has  the standard  observables as well as correlations of kinetic energy, quantum spin-Hall, s-wave pairing and correlations of the generators of SU(2). Note that  the potential and  total  energies  are   defined  as in the  Hamiltonian.   That is   the file Ener_scalJ   corresponds to
  $\langle \hat{H}   \rangle$  with  $\hat{H}$  defined  as  in the  first  equation.

## Limitations
As  it  stands  the  trial wave function is not implemented such that only the finite temperature and not the projective code can be used. Therefore, the code only works for `Projector = .F.` under  `&VAR_Model_Generic`.

```{bibliography}
```