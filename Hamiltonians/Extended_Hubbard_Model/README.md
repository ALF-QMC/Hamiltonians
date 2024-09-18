# Extended  Hubbard  model


## Model definition

The model is defined  as 

$$
\hat{H}  =   \sum_{ b =  \left< i,j \right>  } \sum_{\sigma=1}^{N}  \left(
  \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j}    \hat{c}^{}_{j,\sigma}  +  H.c. \right)   +  \sum_{ b=\left< i,j \right> } \frac{V_{i,j}}{2N} \left( [n_{i} - N/2]  + [ n_{j}  - N/2 ]  \right)^2   + \sum_{i} \frac{U_i}{N} \left( n_i -  \frac{N}{2} \right)^2
$$

The  implementation supports all  standard ALF lattices and  the  hopping as  well as  the
interaction $V$  are  restricted to nearest neighbors.    The parameter file  for  this
specific model reads:

```fortran
&VAR_extended_Hubbard      !! Variables for the Extended Hubbard
ham_T      = 1.d0          ! Hopping parameter
Ham_chem   = 0.d0          ! Chemical potential
Ham_U      = 1.d0          ! Hubbard interaction
Ham_V1     = 0.d0          ! nearest neighbor interaction
ham_T2     = 1.d0          ! For bilayer systems
Ham_U2     = 0.d0          ! For bilayer systems
Ham_V2     = 0.d0          ! For bilayer systems
ham_Tperp  = 0.d0          ! For bilayer systems
Ham_Vperp  = 0.d0          ! For bilayer systems
Ham_Mz     = .T.           ! Logical the choice of the HS field
/
```
In the above Ham_T, Ham_V, Ham_U are the    nearest  neighbor  hopping, nearest  neighbor  interaction  and Hubbard repulsion on the the first layer.   Ham_T2, Ham_V2, Ham_U2   are  the corresponding quantities  on the second layer.  Ham_Tperp, Ham_Vperp   define  the inter-layer couplings.  Finally Ham_chem  is  the chemical potential.
To use  this  Hamiltonian  you  have to specify:
```fortran
&VAR_ham_name
ham_name = "extended_Hubbard"
/
```
in the parameters  file.

In this formulation,  there is  some   double counting  of  the  Hubbard  term.  In particular  expanding  the  square  V-term  gives  the  Hamiltonian:

$$
\hat{H}  =       \sum_{b=\langle  i,j \rangle } \sum_{\sigma=1}^{N}
    \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j}    \hat{c}^{}_{j,\sigma}  +
    \sum_{b=\langle  i,j \rangle }
    \frac{V_{i,j}}{N} \left( [n_{i} - N/2] \right) \left([ n_{j}  - N/2 ]  \right)   +
    \sum_{i}  \frac{U^{eff}_i}{N} \left( n_i -  \frac{N}{2} \right)^2
$$

In the  above

$$
U^{eff}_{i}    = U_i +  \sum_{b = (n,m)} \frac{V_b}{2}  \left( \delta_{i,m} +  \delta_{i,n}\right)
$$

For  the  single  layer  lattices  with uniform $U$ and $V$,  $U^{eff} =   U  + Z V/2$
where  Z=4 (Z=3)   for  the square (honeycomb)  lattice.

If  ```Ham_Mz = .T. ```  the  the code  only works for SU(2) (N=2) spin symmetry and the Haniltonian is written as

$$
\hat{H}  =       \sum_{b=\langle  i,j \rangle } \sum_{\sigma=1}^{N}
    \hat{c}^{\dagger}_{i,\sigma} T^{\sigma}_{i,j}    \hat{c}^{}_{j,\sigma}  +
    \sum_{b=\langle  i,j \rangle }
    \frac{-V_{i,j}}{2N} \sum_{\sigma,\sigma'} \left( [n_{i,\sigma} -1/2]- [n_{i,\sigma'} -1/2]   \right)^2   +
    \sum_{i}  \frac{U^{eff}_i}{N} \left( n_i -  \frac{N}{2} \right)^2
$$

up to a constant, and on each bond we use four fields  to  decouple the $V$-term.  This formulation does better for  sign problem full cases. 

If  `Ham_Mz = .F.`  then the code is based  on   the first equation of the short discription.  We use only one field per bond.  In this case the sign problem is much more severe. 

## Observables 
The code  has  the standard  observables. Note that  the potential and  total  energies  are   defined  as in the  Hamiltonian.   That is   the file Ener_scalJ   corresponds to
  $\langle \hat{H}   \rangle$  with  $\hat{H}$  defined  as  in the  first  (`Ham_Mz = .F.`)  and last (`Ham_Mz = .T.`) equations.

## Limitations
As  it  stands  the  ```Ham_Mz = .T.```  code  is only programmed for  the SU(2) case.  

### Tests

On an $L=6$  square lattice  we  have run the  `Ham_Mz = .T.`   at `Ham_U = 1.5`  and `Ham_V = 0.75` and the `Ham_Mz = .F.`  code  at `Ham_V = 0.75` we  obtain  $-55.023 \pm 0.007$  and $-55.033 \pm 0.001$ for the kinetic energy respectively.  