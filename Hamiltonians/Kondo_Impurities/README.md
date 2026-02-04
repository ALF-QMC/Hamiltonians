# SU(2) Kondo and Anderson impurities.  

**The code can solve the following models**  
**Anderson**

$H = \sum_{i,j} \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j} \hat{c}^{\phantom\dagger}_{j,\sigma} + \sum_{R,R'=1}^{N_{Imp}} t_{R,R'} \hat{d}^{\dagger}_{R,\sigma}\hat{d}^{\phantom\dagger}_{R',\sigma}  + \sum_{R,\sigma} V_R \hat{c}^{\dagger}_{R,\sigma} \hat{d}^{\phantom\dagger}_{R,\sigma}$

**Kondo**

**Name of Hamiltonian: Hamiltonian_Kondo_impurities_smod.F90**

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
- [ ] 

### Usage of the code

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

Imp_t(I,J) Hopping (Anderson) or Heisenberg (Kond) interaction between impurities I,J

Imp_Jz(I,J) Just for Kondo. J_z interaction between I,J impurities

Imp_V(I,R_x,R_y) Hybridization (Anderson) Kondo coupling (Kondo) between impurity I and conduction electron c\_{R_x,R_y}

