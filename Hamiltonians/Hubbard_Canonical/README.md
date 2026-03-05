# Hubbard model in canonical ensemble.

This Hamiltonian extends  the standard implementation of the Hubbard model  with energetically imposed constraints that for calculations in the canonical ensemble  and  fixed z-component of total spin. 

$$\hat{H} = \hat{H}_{\text{Hub}}  +  \lambda_C \left( \hat{N}  -   N_{Particle}\right)^2 + 
    \lambda_S \left( \hat{S}^{z} - S^{z}\right)^2$$

In the above, 
$$\hat{N} =  \sum_{i,n,\delta} \hat{c}^{\dagger}_{i,n,\delta} 
\hat{c}^{}_{i,n,\delta}$$ 
and 
$$\hat{S}^z =  \frac{1}{2}\sum_{i,n,n',\delta} \hat{c}^{\dagger}_{i,n,\delta} 
\sigma^{z}_{n,n'}\hat{c}^{}_{i,n',\delta}$$ 

The  default values  of  $N_{Particle}$ and  $S^z$  are  
$S^z = 0 $ and $N_{Particle} = N * N_{Orb} $  corresponding to  half-filling.  Here $N$ corresponds  to the number of unit cells and $N_{Orb}$  the number of orbitals per  unit cell.

In the present  implementation,  the projection onto the  $S^z =0 $ subspace generates a  mild  sign problem.   For the $4$-site chain at $U/t=4$ and  $\beta t = 2$  the  QMC energy in the $S^z=0$ and $N=4$ sector gives $ -1.8137  \pm  0.0035 $ and the ED  result is $-1.8100613$  so that the present version of the code seems to work. 


Here are  the required   entries for the  parameter file.
```fortran
&VAR_ham_name
ham_name = "Hubbard_Can"
/
```

and 

```fortran
&VAR_Hubbard_Can            !! Variables for the specific model
Mz         = .F.            ! When true, sets the M_z-Hubbard model: Nf=2, demands that
                            ! N_sun is even, HS field couples to the z-component of
                            ! magnetization; otherwise, HS field couples to the density
Continuous = .F.            ! Uses (T: continuous; F: discrete) HS transformation
ham_T      = 1.d0           ! Hopping parameter
ham_chem   = 0.d0           ! Chemical potential
ham_U      = 4.d0           ! Hubbard interaction
ham_Lambda_c = 5.00         ! Projection Charge
ham_Lambda_s = 20.00        ! Projection Spin
ham_T2     = 1.d0           ! For bilayer systems
ham_U2     = 4.d0           ! For bilayer systems
ham_Tperp  = 1.d0           ! For bilayer systems
/
```