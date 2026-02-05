# SU(2) Kondo and Anderson impurities.  

**The code can solve the following models**  
**Anderson**

$H = \sum_{i,j} \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j} \hat{c}^{\phantom\dagger}_{j,\sigma} + \sum_{R,R'=1}^{N_{Imp}} t_{R,R'} \hat{d}^{\dagger}_{R,\sigma}\hat{d}^{\phantom\dagger}_{R',\sigma}  + \sum_{R,\sigma} V_R \hat{c}^{\dagger}_{R,\sigma} \hat{d}^{\phantom\dagger}_{R,\sigma}  + \text{H.c.}  + \frac{U}{2} \sum_{R}\left( \hat{n}^{d}_R - 1 \right)^2$

**Kondo**

$H = \sum_{i,j} \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j} \hat{c}^{\phantom\dagger}_{j,\sigma} + \sum_{R,R'=1}^{N_{Imp}} J_{R,R'} \vec{\hat{S}}_{R} \cdot \vec{\hat{S}}_{R'}  +  
\sum_{R,R'=1}^{N_{Imp}}  J^{z}_{R,R'} \hat{S}^{z}_{R}  \hat{S}^{z}_{R'}  + \sum_{R} J^{K}_{R} \vec{\hat{S}}_{R}  \cdot   \vec{\hat{c}}^{\dagger}_{R} \frac{\vec{\sigma}}{2} \vec{\hat{c}}^{\phantom\dagger}_{R} $ 

In the above the hopping is on a square lattice  with the possibility of an altermagnetic term  as  discussed in [this](https://arxiv.org/abs/2601.07138) paper.  

**Name of Hamiltonian: Hamiltonian_Kondo_impurities_smod.F90**

Here are the of things that are implemented: 

- [x] Symmetric Trotter
- [x] Anderson as well as Kondo impurities
- [x] Specify impurity directly in the parameter file
- [x] Spin correlations between the impurities
- [x] Local composite fermion (Kondo) or f-(Anderson) Green function.
- [x] We again need two lattices. For the impurities the lattice has one unit cell and the number of orbitals correspond to the number of impurities. For the conduction electrons we can use the square lattice.
- [x] Easy axis anisotropy (Sz Sz) correlations for local spins. Just for Kondo.
- [x] Particle-hole symmetry is automatically detected and if true, flavor symmetry is used.
- [x] Added possibility of alter-magnetic  or spin-nematic band.  $d_{xy}$ and $d_{x^2-y^2}$
- [ ] Magnetic field. For Kondo as well as for Anderson.  This really not very hard and could be implemented quickly
- [ ] Projector

Please do not hesitate to contact us if you need more features.

### Usage of the code

## Input

Here are examples for parameters.  


Hopping and Hubbard U
```
&VAR_Kondo_impurities    !! Variables for the Kondo impurity  code
ham_T     = 1.d0            ! Hopping parameter !  Conduction.
ham_chem  = 0.d0            ! Chemical potential
Ham_T_alterm = 0.40         ! Alter-magnetic hopping just conduction.
d_wave_rep   ="DXY"         ! DXY or DX2Y2.Representation of d-wave alter-magnetic hopping. 
ham_U     = 2.d0            ! Hubbard
Ham_Imp_Kind = "Anderson"   ! Impurity kind "Anderson" or   "Kondo" 
Ham_N_imp = 6               ! # Spin 1/2  impurities 
/
```
Position and interaction of impurities
```
&VAR_impurities    !! Variables for the Kondo impurity  code
Imp_t(1,2)    =  -1.d0 
Imp_t(2,3)    =  -1.d0
Imp_t(3,4)    =  -1.d0
Imp_t(4,5)    =  -1.d0
Imp_t(5,6)    =  -1.d0
Imp_t(6,1)    =  -1.d0
Imp_Jz(1,2)  =  2.d0
Imp_Jz(1,3)  =  2.d0
Imp_Jz(2,3)  =  2.d0
Imp_V(1,1,1)  =  1.d0
Imp_V(2,1,2)  =  1.d0
Imp_V(3,1,3)  =  1.d0
/
```
The input parameters have different meanings for different model choice.

**Anderson**

Imp_V(n,R_x,R_y) Hybridization between impurity n and conduction electron 
$\hat{c}_{R_x,R_y,\sigma}$

Imp_t(n,m) Hopping between impurities  n and m

**Kondo** 

Imp_V(n,R_x,R_y) Kondo coupling between impurity n and conduction electron 
$\hat{c}_{R_x,R_y,\sigma}$

Imp_t(n,m) Heisenberg interaction between impurities n,m

Imp_Jz(n,m) J_z interaction between n,m impurities


## Output

As it stands the program computes the local site dependent [composite fermion  Green function](https://arxiv.org/abs/2107.10272)  and the  spin-spin correlations.  

The position and numbering of the impurities is given by the Imp_V(n,Rx,Ry)   input data.

**SU(2) code**   (N_SUN = 2, N_FL = 1)

GreenPsi_n_Rx_Ry_tau   contains  the local composite fermion Green function  $\langle \Psi_{R}(\tau)\Psi^{\dagger}_{R}(\tau) \rangle $ 

SpinZ_tau  contains  the  dynamical spin-spin correlations  $S(n,m,\tau) = \langle \hat{S}^z_n(\tau) \hat{S}^z_m(0) \rangle$  

**U(1) code** (N_SUN = 1, N_FL = 2)   

In the presence of the altermagnetic term,   SU(2) symmetry is broken down to U(1)  and  the files mentioned above acquire a spin index, up,  down.  Aside from this, the structure is the same. 

## Analysis 
To analyze  the  data  run the  analysis code  $ALF_DIR/Analysis/ana.out *  (or ana_hdf5 if you are opting for HDF5 formatted output.)
