!  Copyright (C) 2016 - 2018 The ALF project
! 
!     The ALF project is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     The ALF project is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with ALF.  If not, see http://www.gnu.org/licenses/.
!
!     Under Section 7 of GPL version 3 we require you to fulfill the following additional terms:
!
!     - It is our hope that this program makes a contribution to the scientific community. Being
!       part of that community we feel that it is reasonable to require you to give an attribution
!       back to the original authors if you have benefitted from this program.
!       Guidelines for a proper citation can be found on the project's homepage
!       http://alf.physik.uni-wuerzburg.de .
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version.

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!> This module provides a set of predefined lattices, hoppings, interactions,
!> trial wave functions as well as observables.
!>       
!
!--------------------------------------------------------------------

    Module Predefined_structures
      use iso_fortran_env, only: output_unit, error_unit
      
      Use Lattices_v3
      Use Operator_mod
      Use WaveFunction_mod
      Use MyMats
      Implicit none

      private

      public :: Predefined_checkerboard, checkeboard_add_subfam, &
                predefined_trialwavefunction, Predefined_Latt
      
      
    contains

!--------------------------------------------------------------------
!> @author
!> ALF-project
!
!> @brief
!> Definition of  a set of lattices: Square,  Honeycomb, Pi_Flux
!>
!> @param [in]  Latttice_type
!>\verbatim
!> Character(64)
!> Can take the values
!> Square,  Honeycomb, Pi_Flux
!> \endverbatim
!> @param [in]  L1, L2
!>\verbatim
!>    Integer
!>    Size of the lattice in units of the lattice constants
!>\endverbatim
!> @param [out]  Latt_unit
!>\verbatim
!>    Type (Unit_cell)
!>    The unit cell. Contains Norb, N_coord and positions of orbitals.
!>\endverbatim
!> @param [out]  Ndim
!>\verbatim
!>    Integer
!>    Number of orbitals
!>\endverbatim
!> @param [out]  List, Invlist
!>\verbatim
!>    Integer(:,:)
!>    List(I=1.. Ndim,1)    =   Unit cell of site I
!>    List(I=1.. Ndim,2)    =   Orbital index  of site I
!>    Invlist(1..Unit_cell,1..Orbital) = site I
!>\endverbatim
!> @param [out]  Latt
!>\verbatim
!>    Type(Lattice)
!>    Sets the lattice
!>\endverbatim
!> @param [out]  Latt_unit
!>\verbatim
!>    Type(Unit_cell)
!>    Sets the lattice
!>\endverbatim
!>
!-------------------------------------------------------------------
      Subroutine Predefined_Latt(Lattice_type, L1, L2, &
            Ndim, List, Invlist, Latt, Latt_Unit, List_orb, Invlist_orb)

        Implicit none

        !Set the lattice
        Character (len=64), Intent(IN)                     :: Lattice_type
        Integer, Intent(IN)                                :: L1,L2
        Integer, Intent(OUT)                               :: Ndim
        Integer, Intent(OUT), Dimension(:,:), allocatable  :: List, Invlist
        Type(Unit_cell), Intent(Out)                       :: Latt_Unit
        Type(Lattice), Intent(Out)                         :: Latt
        Integer, Intent(inout), Dimension(:,:), allocatable  :: List_orb, Invlist_orb
        Real (Kind=Kind(0.d0))  :: A1_p(2), a2_p(2), L1_p(2), L2_p(2)
        Integer :: I, nc, no,n, no1, no2

        select case (Lattice_type)
        case("Square")
          If ( L2==1 .and. L1 > 1 ) then
            Latt_Unit%N_coord   = 1
          elseif (L2 >1 .and. L1 > 1) then
            Latt_Unit%N_coord   = 2
          else
            Write(error_unit,*) 'For one-dimensional lattices set L2=1.'
            Write(error_unit,*) 'You can also use use n_leg_ladder with n=1'
            error stop 1
          endif
          Latt_Unit%Norb      = 1
          Allocate (Latt_unit%Orb_pos_p(1,2))
          allocate (List_orb(0, 0), Invlist_orb(0, 0))
          Latt_Unit%Orb_pos_p(1,:) = 0.d0
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =  dble(L2)*a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt )
        case("N_leg_ladder")
          a1_p(1) =  1.0  ; a1_p(2) =  0.d0
          a2_p(1) =  0.0  ; a2_p(2) =  1.d0
          L1_p    =  dble(L1)*a1_p
          L2_p    =           a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

          Latt_Unit%Norb     = L2
          Latt_Unit%N_coord  = 1
          Allocate (Latt_unit%Orb_pos_p(L2,2))
          allocate (List_orb(0, 0), Invlist_orb(0, 0))
          do no = 1,L2
            Latt_Unit%Orb_pos_p(no,1) = 0.d0
            Latt_Unit%Orb_pos_p(no,2) = real(no-1,kind(0.d0))
          enddo
        case("Open_square")
          a1_p(1) = real(L1, kind(0.d0)); a1_p(2) = 0.d0
          a2_p(1) = 0.d0; a2_p(2) = real(L2, kind(0.d0))
          L1_p = a1_p
          L2_p = a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p, a2_p, Latt )

          Latt_Unit%Norb     = L1*L2
          Latt_Unit%N_coord  = 0
          Allocate (Latt_unit%Orb_pos_p(Latt_Unit%Norb,2))
          allocate (List_orb(Latt_Unit%Norb, 2), Invlist_orb(L1, L2))
          no = 0
          do no1 = 1, L1
            do no2 = 1, L2
              no = no + 1
              Latt_Unit%Orb_pos_p(no,1) = real(no1-1,kind(0.d0))
              Latt_Unit%Orb_pos_p(no,2) = real(no2-1,kind(0.d0))
              List_orb(no,1) = no1
              List_orb(no,2) = no2
              Invlist_orb(no1, no2) = no
            enddo
          enddo
        case default
          Write(error_unit,*) "Predefined_Latt: Lattice not yet implemented!"
          error stop 1
        end select

        Ndim = Latt%N*Latt_Unit%Norb
        Allocate (List(Ndim,2), Invlist(Latt%N,Latt_Unit%Norb))
        nc = 0
        Do I = 1,Latt%N
          Do no = 1,Latt_Unit%Norb
              ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
              nc = nc + 1
              List(nc,1) = I
              List(nc,2) = no
              Invlist(I,no) = nc
          Enddo
        Enddo

      end Subroutine Predefined_Latt

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief
      !> Adds a second index to the orbital structure of a lattice
!>
!> @param [in]  Latt
!>\verbatim 
!>    Type(Lattice)
!>    The bravais lattice
!>\endverbatim 
!> @param [in]  Norb
!>\verbatim 
!>    Integer  
!>    Number of orbitals per unit cell due to lattice geometry
!>\endverbatim 
!> @param [in]  Norb2
!>\verbatim 
!>    Integer  
!>    Number of sub-orbitals per lattice orbital 
!>\endverbatim 
!> @param [out]  Ndim
!>\verbatim 
!>    Integer
!>    Number of orbitals      
!>\endverbatim 
!> @param [out]  List, Invlist
!>\verbatim 
!>    Integer(:,:)
!>    List(I=1.. Ndim,1)    =   Unit cell of site I    
!>    List(I=1.. Ndim,2)    =   First orbital index  of site I  
!>    List(I=1.. Ndim,3)    =   Second orbital index  of site I    
!>    Invlist(1..Unit_cell,1..Orbital1,,1..Orbital2) = site I    
!>\endverbatim 
!>
!-------------------------------------------------------------------
      Subroutine Latt_add_orbitals(N_unit, N_orb, N_orb2, List, Invlist)

        Implicit none

        !Set the lattice
        Integer, Intent(IN)                :: N_unit, N_orb, N_orb2
        Integer, Intent(OUT), allocatable  :: List(:,:), Invlist(:,:,:)
        
        Integer :: I, nc, no, no2, Ndim

        Ndim = N_unit*N_orb*N_orb2
        Allocate (List(Ndim,3), Invlist(N_unit,N_orb,N_orb2))
        nc = 0
        Do I = 1,N_unit
           Do no = 1,N_orb
              ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
              Do no2 = 1,N_orb2
                 nc = nc + 1
                 List(nc,1) = I
                 List(nc,2) = no
                 List(nc,3) = no2
                 Invlist(I,no,no2) = nc 
              Enddo
           Enddo
        Enddo
      end Subroutine Latt_add_orbitals

      
      Subroutine Hopping_add_orbitals(OP_old, OP_new, N_sub)
         Implicit none
         Type(Operator), Intent(In)               :: OP_old(:,:)
         Type(Operator), Intent(Out), allocatable :: OP_new(:,:)
         Integer, Intent(IN)                      :: N_sub
         
         Integer :: I, J, nI, nJ, I2, J2, Ndim_old, Ndim_new, n
         Complex (Kind=Kind(0.d0)) :: chem, t
      
         Ndim_old = size(OP_old(1,1)%O,1)
         
         allocate(Op_new(1,1))
         Ndim_new = Ndim_old*N_sub
         
         Call Op_make(Op_new(1,1),Ndim_new)
         
         do I=1,Ndim_old
            do J=1,Ndim_old
               t = OP_old(1,1)%O(I,J)
               do n=1, N_sub
                  I2 = (I-1)*N_sub + n
                  J2 = (J-1)*N_sub + n
                  OP_new(1,1)%O(I2,J2) = t
               enddo
            enddo
         enddo
         
         Do I = 1,Ndim_old
            do nI = 1, N_sub
               I2 = (I-1)*N_sub + nI
               Op_new(1,1)%P(I2) = I2
            enddo
         Enddo
         Op_new(1,1)%g     = Op_old(1,1)%g
         Op_new(1,1)%alpha = Op_old(1,1)%alpha
         Call Op_set(Op_new(1,1))

      end Subroutine Hopping_add_orbitals

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!>
!> @brief 
!>    Hopping, with or without checkerboard  
!>    Per flavor, the  hopping is given by
!>    \f[  e^{ - \Delta \tau  H_t  }   = \prod_{n=1}^{N} e^{ - \Delta \tau_n  H_t(n) }   \f]
!>    If  _Symm_ is set to true and if  _Checkeborad_ is on, then  one will carry out a
!>    symmetric decomposition so as to preserve  the hermitian properties of the hopping.
!>    Thereby   OP_T has dimension OP_T(N,N_FL)
!> @param [in]  Latttice_type
!>    Character(64)
!>\verbatim 
!>     Square,  Honeycomb, Pi_Flux 
!>\endverbatim 
!> @param [in]  Latt_unit
!>    Type(Unit_cell)
!> \verbatim
!>     Contains number of orbitals per unit cell and positions, as well as coordination number
!> \endverbatim
!> @param [in]  Ndim
!>    Integer
!> \verbatim
!>     Number of orbitals      
!> \endverbatim
!> @param [in]  List, Invlist
!>    Integer(:,:)
!> \verbatim
!>      List(I=1.. Ndim,1)    =   Unit cell of site I    
!>      List(I=1.. Ndim,2)    =   Orbital index  of site I    
!>      Invlist(Unit_cell,Orbital) = site I    
!> \endverbatim
!> @param [in]    Latt
!>    Type(Lattice)
!> \verbatim
!>      The Lattice
!> \endverbatim
!> @param [in]  Dtau
!>    Real
!> \verbatim
!>      Imaginary time step
!> \endverbatim
!> @param [in]  Ham_T
!>    Real
!> \verbatim
!>      Hopping matrix element
!> \endverbatim
!> @param [in]  Ham_Chem
!>    Real
!> \verbatim
!>      Chemical potential
!> \endverbatim
!> @param [in]  XB_X, YB_Y
!>    Real
!> \verbatim
!>      X, Y  Boundary conditions
!> \endverbatim
!> @param [in]  Phi_X, Phi_Y 
!>    Real
!> \verbatim
!>      X, Y  Fluxes
!> \endverbatim
!> @param [in]  N_FL
!>    Integer
!> \verbatim
!>      Flavors
!> \endverbatim
!> @param [in]  Checkerboard
!>    Logical
!> \verbatim
!>      Allows for checkerboard decomposition
!> \endverbatim
!> @param [in]  Symm
!>    Logical
!> \verbatim
!>      Allows for symmetric checkerboard decomposition
!> \endverbatim
!> @param [out]  OP_T 
!>    Type(operator)(N,N_FL)
!> \verbatim
!>      Hopping
!> \endverbatim
!> @param [in]  Dimer 
!>    Real, Optional.  Modulation of hopping that breaks lattice symmetries so as to generate a unique
!>    ground state for the half-filled case.  This option is  effective only  if the checkerboard
!>    decomposition is not used. It is presently implemented for the square and one-dimensional lattices.
!> \verbatim
!>      Hopping
!> \endverbatim
!>       
!------------------------------------------------------------------
      Subroutine Predefined_Hopping0(Lattice_type, Ndim, List,Invlist,Latt,  Latt_unit,  &
         &                        L1, L2, List_orb, Invlist_orb, &
         &                        Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
         &                        N_FL, OP_T )

        Implicit none

        Character (len=64), Intent(IN)               :: Lattice_type
        Integer, Intent(IN)                          :: Ndim, N_FL, L1, L2
        Integer, Intent(inout), Dimension(:,:)          :: List, Invlist, List_orb, Invlist_orb
        Type(Lattice),  Intent(in)                   :: Latt
        Type(Unit_cell),Intent(in)                   :: Latt_unit
        Real (Kind=Kind(0.d0)), Intent(In)           :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y
        
        Type(Operator), Intent(Out),  dimension(:,:), allocatable  :: Op_T 


        !Local
        Integer :: I, I1, J1, I2, n, Ncheck,nc, nc1, no, N_Fam, L_FAM, I1x, I1y, n1, I2x, I2y
        Complex (Kind=Kind(0.d0)) :: ZX, ZY
        Real    (Kind=Kind(0.d0)) :: del_p(2), X, g


         allocate(Op_T(1,N_FL))
         do n = 1,N_FL
            Call Op_make(Op_T(1,n),Ndim)
            nc = 1
            Select case (Lattice_type)
            Case ("Square")
              ZX  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/Xnorm(Latt%L1_p), kind=kind(0.d0) ) )
              ZY  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_Y/Xnorm(Latt%L2_p), kind=kind(0.d0) ) )
              DO I = 1, Latt%N
                I1 = Latt%nnlist(I,1,0)
                I2 = Latt%nnlist(I,0,1)
                If ( Latt%list(I,1) == 0 ) then
                  Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*ZX
                  Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(ZX)
                else
                  Op_T(nc,n)%O(I,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))*ZX
                  Op_T(nc,n)%O(I1,I) = cmplx(-Ham_T, 0.d0, kind(0.D0))*conjg(ZX)
                endif
                If ( Latt%list(I,2) == 0 ) then
                  Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T*XB_Y,    0.d0, kind(0.D0))*ZY
                  Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T*XB_Y,    0.d0, kind(0.D0))*conjg(ZY)
                else
                  Op_T(nc,n)%O(I,I2) = cmplx(-Ham_T     ,    0.d0, kind(0.D0))*ZY
                  Op_T(nc,n)%O(I2,I) = cmplx(-Ham_T     ,    0.d0, kind(0.D0))*conjg(ZY)
                endif
                Op_T(nc,n)%O(I ,I) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
              enddo
            Case ("N_leg_ladder")
              ZX  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/Xnorm(Latt%L1_p), kind=kind(0.d0) ) )
              do I1x = 1, Latt%N
                do I1y = 1, Latt_Unit%Norb
                  I1 = invlist(I1x, I1y)
                  Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  do nc1 = 1, 2
                    I2x = I1x
                    I2y = I1y
                    if (nc1 == 1 ) I2x = I1x+1
                    if (nc1 == 2 ) I2y = I1y+1
                    if (I2y > Latt_Unit%Norb) cycle
                    if (I2x > Latt%N) then
                      I2x = 1
                      I2 = invlist(I2x, I2y)
                      Op_T(nc,n)%O(I1,I2) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*ZX
                      Op_T(nc,n)%O(I2,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(ZX)
                    else
                      I2 = invlist(I2x, I2y)
                      Op_T(nc,n)%O(I1,I2) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                    endif
                  enddo
                enddo
              enddo
            Case ("Open_square")
              ZX  =  exp( cmplx(0.d0, 2.d0 * acos(-1.d0)*Phi_X/Xnorm(Latt%L1_p), kind=kind(0.d0) ) )
              do I1x = 1, L1
                do I1y = 1, L2
                  I1 = Invlist_orb(I1x, I1y)
                  Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                  do nc1 = 1, 2
                    I2x = I1x
                    I2y = I1y
                    if (nc1 == 1 ) I2x = I1x+1
                    if (nc1 == 2 ) I2y = I1y+1
                    ! if (I2x > L1 .or. I2y > L2) cycle
                    ! I2 = Invlist_orb(I2x, I2y)
                    ! Op_T(nc,n)%O(I1,I2) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                    ! Op_T(nc,n)%O(I2,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                    if (I2y > L2) cycle
                    if (I2x > L1) then
                      I2x = 1
                      I2 = Invlist_orb(I2x, I2y)
                      Op_T(nc,n)%O(I1,I2) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*ZX
                      Op_T(nc,n)%O(I2,I1) = cmplx(-Ham_T*XB_X, 0.d0, kind(0.D0))*conjg(ZX)
                    else
                      I2 = Invlist_orb(I2x, I2y)
                      Op_T(nc,n)%O(I1,I2) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                      Op_T(nc,n)%O(I2,I1) = cmplx(-Ham_T, 0.d0, kind(0.D0))
                    endif
                  enddo
                enddo
              enddo
            Case ("Honeycomb")
              X = 2.d0 * acos(-1.d0)*Phi_X / ( Xnorm(Latt%L1_p) * (Xnorm(Latt%a1_p)**2)  )
              DO I = 1, Latt%N
                do no = 1,Latt_unit%Norb
                  I1 = Invlist(I,no)
                  Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                enddo
                I1 = Invlist(I,1)
                J1 = I1
                Do nc1 = 1,Latt_unit%N_coord
                  select case (nc1)
                  case (1)
                    J1 = invlist(I,2)
                    del_p(:)  =  Latt_unit%Orb_pos_p(2,:) 
                  case (2)
                    J1 = invlist(Latt%nnlist(I,1,-1),2)
                    del_p(:)   =  Latt%a1_p(:) - Latt%a2_p(:)  + Latt_unit%Orb_pos_p(2,:)
                  case (3)
                    J1 = invlist(Latt%nnlist(I,0,-1),2) 
                    del_p(:)   =  - Latt%a2_p(:) +  Latt_unit%Orb_pos_p(2,:) 
                  case default
                    Write(error_unit,*) ' Error in  Predefined_Hopping0 '  
                    error stop
                  end select
                  ZX = exp( cmplx(0.d0, X*Iscalar(Latt%a1_p,del_p), kind(0.D0) ) )
                  Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * ZX
                  Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0)) * CONJG(ZX) 
                Enddo
              Enddo
            case("Pi_Flux")
              DO I = 1, Latt%N
                do no = 1,Latt_unit%Norb
                  I1 = Invlist(I,no)
                  Op_T(nc,n)%O(I1 ,I1) = cmplx(-Ham_chem, 0.d0, kind(0.D0))
                enddo
                I1 = Invlist(I,1)
                J1 = I1
                Do nc1 = 1,Latt_unit%N_coord
                  select case (nc1)
                  case (1)
                    J1 = invlist(I,2) 
                  case (2)
                    J1 = invlist(Latt%nnlist(I,0, 1),2) 
                  case (3)
                    J1 = invlist(Latt%nnlist(I,-1,1),2) 
                  case (4)
                    J1 = invlist(Latt%nnlist(I,-1,0),2) 
                  case default
                    Write(error_unit,*) ' Error in Predefined_Hopping0 '  
                    error stop
                  end select
                  if (nc1 == 1 ) then
                    Op_T(nc,n)%O(I1,J1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                    Op_T(nc,n)%O(J1,I1) = cmplx( Ham_T,    0.d0, kind(0.D0))
                  Else
                    Op_T(nc,n)%O(I1,J1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                    Op_T(nc,n)%O(J1,I1) = cmplx(-Ham_T,    0.d0, kind(0.D0))
                  endif
                Enddo
              Enddo
            case default 
              Write(error_unit,*) "Predefined_Hopping0: Lattice '", Lattice_type, "' not yet implemented!"
              error stop
            end select
            Do I = 1,Ndim
              Op_T(nc,n)%P(i) = i 
            Enddo
            if ( abs(Ham_T) < 1.E-6  .and.  abs(Ham_chem) < 1.E-6 ) then 
              Op_T(nc,n)%g = 0.d0
            else
              Op_T(nc,n)%g = -Dtau
            endif
            Op_T(nc,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
            Call Op_set(Op_T(nc,n))
            !Do I = 1,Size(Op_T(nc,n)%E,1)
            !   Write(6,*) Op_T(nc,n)%E(I)
            !Enddo
         Enddo
      
    End Subroutine Predefined_Hopping0


!--------------------------------------------------------------------
!> @author 
      !> ALF-project
      !>
      !> @brief 
      !>    Sets the trial wave function corresponding to the solution of the non-interacting
      !>    tight binding Hamiltonian on the given lattice. Twisted boundary conditions (Phi_X=0.01)
      !>    are implemented so as to generate a non-degenerate trial wave functions. 
      !> @param [in]  Lattice_type
      !>    Character(64)
      !> \verbatim
      !>    Square,  Honeycomb, Pi_Flux
      !> \endverbatim
      !> @param [in]  Latt_unit
      !>    Type(Unit_cell)
      !> \verbatim
      !>     Contains number of orbitals per unit cell and positions, as well as coordination number
      !> \endverbatim
      !> @param [in]  Ndim
      !>    Integer
      !> \verbatim
      !>     Number of orbitals
      !> \endverbatim
      !> @param [in]  List, Invlist
      !>    Integer(:,:)
      !> \verbatim
      !>    List(I=1.. Ndim,1)    =   Unit cell of site I    
      !>    List(I=1.. Ndim,2)    =   Orbital index  of site I    
      !>    Invlist(Unit_cell,Orbital) = site I    
      !> \endverbatim
      !> @param [in]    Latt
      !>    Type(Lattice)
      !> \verbatim
      !>    The Lattice
      !> \endverbatim
      !> @param [in]  N_part
      !>    Integer
      !> \verbatim
      !>    Particle number for each flavor
      !> \endverbatim
      !> @param [in]  N_FL
      !>    Integer
      !> \verbatim
      !>    Flavor
      !> \endverbatim
      !> @param [out]  WF_L, WF_R
      !>    Type(Wavefunction)(N_FL)
      !> \verbatim
      !>    Wavefunction
      !>    Also sets the degeneracy:  E(N_part + 1) - E(N_part). Energy eigenvalues are ordered in ascending order.
      !> \endverbatim
      !>       
!------------------------------------------------------------------       
      Subroutine Predefined_TrialWaveFunction(Lattice_type, L1, L2, N_part, N_FL, WF_L, WF_R, N_sub)

        Implicit none
        Character (len=64), Intent(IN)                     :: Lattice_type
        Integer, Intent(IN)                                :: L1, L2, N_part
        Type(WaveFunction), Intent(out), Dimension(:), allocatable :: WF_L, WF_R
        Integer, Intent(in), optional :: N_sub
        
        Integer              :: Ndim, N_FL
        Integer, allocatable :: List(:,:), Invlist(:,:), List_orb(:,:), Invlist_orb(:,:)
        Type(Lattice)        :: Latt
        Type(Unit_cell)      :: Latt_Unit

        
        Type(Operator),  dimension(:,:), allocatable  :: OP_tmp, OP_tmp0
        Real (Kind=Kind(0.d0))                        :: Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, Dimer

        Integer :: N, nf, I1, I2
        
        Dtau     = 1.d0
        Ham_T    = 1.d0
        Ham_Chem = 0.d0
        XB_X     = 1.d0
        XB_Y     = 1.d0
        Phi_X    = 0.01
        Phi_Y    = 0.d0
        Dimer    = 0.d0
        
        call Predefined_Latt(Lattice_type, L1, L2, &
              Ndim, List, Invlist, Latt, Latt_Unit, List_orb, Invlist_orb)
        if ( present(N_sub) ) then
          Call Predefined_Hopping0(Lattice_type, Ndim, List, Invlist, Latt, Latt_Unit, &
            &                    L1, L2, List_orb, Invlist_orb, &
            &                    Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
            &                    N_FL, OP_tmp0 )
          call Hopping_add_orbitals(OP_tmp0, OP_tmp, N_sub)
          Ndim = Ndim*N_sub
        else
          Call Predefined_Hopping0(Lattice_type ,Ndim, List,Invlist,Latt, Latt_Unit, &
              &                    L1, L2, List_orb, Invlist_orb, &
              &                    Dtau, Ham_T, Ham_Chem, XB_X, XB_Y, Phi_X, Phi_Y, &
              &                    N_FL, OP_tmp )
        endif
        
        Allocate(WF_L(N_FL),WF_R(N_FL))
        do n=1,N_FL
          Call WF_alloc(WF_L(n),Ndim,N_part)
          Call WF_alloc(WF_R(n),Ndim,N_part)
        enddo

        Do nf = 1,N_FL
          Call Diag(Op_tmp(1,nf)%O,Op_tmp(1,nf)%U,Op_tmp(1,nf)%E)
          do I2=1,N_part
              do I1=1,Ndim
                WF_L(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
                WF_R(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
              enddo
          enddo
          WF_L(nf)%Degen = Op_tmp(1,1)%E(N_part+1) - Op_tmp(1,1)%E(N_part)
          WF_R(nf)%Degen = Op_tmp(1,1)%E(N_part+1) - Op_tmp(1,1)%E(N_part)
        enddo
         
        Do nf = 1,N_FL
          Call Op_clear(OP_tmp(1,nf),Ndim)
        enddo
        Deallocate (OP_tmp)
        if (allocated(OP_tmp0)) then
          Do nf = 1,N_FL
            Call Op_clear(OP_tmp0(1,nf),Ndim)
          enddo
          Deallocate (OP_tmp0)
        endif
        Deallocate (List, Invlist)
         
      end Subroutine Predefined_TrialWaveFunction

      
      !> checker(3, L_Fam, N_Fam)
      !> first index 1,2: site indices of bonds
      !> first index 3: type of bond
      !>             default: 0
      !>           x-boundary => checker(3, .., ..) += 1
      !>           y-boundary => checker(3, .., ..) += 2
      !>           dimer-bond => checker(3, .., ..) += 4
      Subroutine Predefined_checkerboard(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
                                    &  L1, L2, List_orb, Invlist_orb, &
                                    &  checker, L_Fam)

        Implicit none

        Character (len=64), Intent(IN)           :: Lattice_type
        Integer, Intent(IN)                      :: Ndim
        Integer, Intent(IN), Dimension(:,:)      :: List, Invlist
        Type(Lattice), Intent(in)                :: Latt
        Type(Unit_cell), Intent(in)              :: Latt_unit
        Integer, Intent(IN)                      :: L1, L2
        Integer, Intent(inout), Dimension(:,:)      :: List_orb, Invlist_orb
        
        Integer, Intent(Out), allocatable :: checker(:,:,:)
        integer, Intent(Out), allocatable :: L_Fam(:)

        !Local
        Integer :: L_Fam_max, N_Fam
        Integer :: I, I1, I2, nc, nc1, nf, I1x, I1y, I2x, I2y
        
        Integer :: N_x_bound, N_y_bound
        
        select case (Lattice_type)
        case ("Square")
          N_Fam = 4             !  Number of Families = 4
          L_Fam_max = Latt%N/2  !  Length of Families = LQ/4
          allocate( L_Fam(N_Fam), checker(2, L_Fam_max, N_Fam) )
          checker(:,:,:) = 0
          nc = 0
          do nc1 = 1,N_Fam
            nf = 0
            L_Fam(nc1) = L_Fam_max
            do I = 1,Latt%N
              if ( mod(Latt%List(I,1) + Latt%List(I,2),2) == 0 ) then
                !if( nf == L_Fam(nc1) ) nf = 0
                nc = nc + 1
                nf = nf + 1
                I1 = I
                if (nc1 == 1 ) I2 = latt%nnlist(I1, 1, 0) 
                if (nc1 == 2 ) I2 = latt%nnlist(I1, 0, 1)
                if (nc1 == 3 ) I2 = latt%nnlist(I1,-1, 0)
                if (nc1 == 4 ) I2 = latt%nnlist(I1, 0,-1)
                checker(1, nf, nc1) = I1
                checker(2, nf, nc1) = I2
              endif
            enddo
          enddo
        case ("N_leg_ladder")
          N_Fam = 4                            !  Number of Families = 4
          L_Fam_max = Latt%N*Latt_Unit%Norb/2  !  Length of Families = LQ/4
          allocate( L_Fam(N_Fam), checker(2, L_Fam_max, N_Fam) )
          checker(:,:,:) = 0
          nc = 0
          do nc1 = 1,N_Fam
            nf = 0
            L_Fam(nc1) = L_Fam_max
            do I1x = 1, Latt%N
              do I1y = 1, Latt_Unit%Norb
                if ( mod(I1x + I1y, 2) == 0 ) then
                  if( nf == L_Fam(nc1) ) nf = 0
                  nc = nc + 1
                  nf = nf + 1
                  I2x = I1x
                  I2y = I1y
                  if (nc1 == 1 ) I2x = I1x+1
                  if (nc1 == 2 ) I2y = I1y+1
                  if (nc1 == 3 ) I2x = I1x-1
                  if (nc1 == 4 ) I2y = I1y-1
                  if (I2y < 1 .or. I2y > Latt_Unit%Norb) then
                    nc = nc - 1
                    nf = nf - 1
                    L_Fam(nc1) = L_Fam(nc1) - 1
                  else
                    if (I2x < 1) I2x = Latt%N
                    if (I2x > Latt%N) I2x = 1
                    checker(1, nf, nc1) = invlist(I1x, I1y)
                    checker(2, nf, nc1) = invlist(I2x, I2y)
                  endif
                endif
              enddo
            enddo
          enddo
        case ("Open_square")
          N_Fam = 4                            !  Number of Families = 4
          L_Fam_max = L1*L2/2  !  Length of Families = LQ/4
          allocate( L_Fam(N_Fam), checker(2, L_Fam_max, N_Fam) )
          checker(:,:,:) = 0
          nc = 0
          do nc1 = 1,N_Fam
            nf = 0
            L_Fam(nc1) = L_Fam_max
            do I1x = 1, L1
              do I1y = 1, L2
                if ( mod(I1x + I1y, 2) == 0 ) then
                  if( nf == L_Fam(nc1) ) nf = 0
                  nc = nc + 1
                  nf = nf + 1
                  I2x = I1x
                  I2y = I1y
                  if (nc1 == 1 ) I2x = I1x+1
                  if (nc1 == 2 ) I2y = I1y+1
                  if (nc1 == 3 ) I2x = I1x-1
                  if (nc1 == 4 ) I2y = I1y-1
                  if (I2x < 1 .or. I2x > L1 .or. I2y < 1 .or. I2y > L2) then
                    nc = nc - 1
                    nf = nf - 1
                    L_Fam(nc1) = L_Fam(nc1) - 1
                  else
                    checker(1, nf, nc1) = Invlist_orb(I1x, I1y)
                    checker(2, nf, nc1) = Invlist_orb(I2x, I2y)
                  endif
                endif
              enddo
            enddo
          enddo
        case ("Honeycomb")
          N_Fam = 3           !  Number of Families = 3
          L_Fam_max = Latt%N  !  Length of Families
          allocate( L_Fam(N_Fam), checker(2, L_Fam_max, N_Fam) )
          nc = 0
          do nc1 = 1, N_Fam
            L_Fam(nc1) = L_Fam_max
            do I = 1, L_Fam_max
              nc = nc + 1
              I1 = invlist(I,1)
              if (nc1 == 1 ) I2 = invlist(I,2)
              if (nc1 == 2 ) I2 = invlist(latt%nnlist(I,1,-1),2)
              if (nc1 == 3 ) I2 = invlist(latt%nnlist(I,0,-1),2)
              checker(1, I, nc1) = I1
              checker(2, I, nc1) = I2
            enddo
          enddo
        case ("One_dimensional")
          N_Fam = 2              !  Number of Families = 2
          L_Fam_max = Latt%N /2  !  Length of Families
          allocate( L_Fam(N_Fam), checker(2, L_Fam_max, N_Fam) )
          nc = 0
          do nc1 = 1, N_Fam
            L_Fam(nc1) = L_Fam_max
            do I = 1, L_Fam_max
              nc = nc + 1
              I1 = I*2
              if (nc1 == 1 ) I2 = latt%nnlist(I1,1 ,0)
              if (nc1 == 2 ) I2 = latt%nnlist(I1,-1,0)
              checker(1, I, nc1) = I1
              checker(2, I, nc1) = I2
            enddo
          enddo
        case default
          Write(error_unit,*) "Predefined_checkerboard is not implemented for this lattice", Lattice_type
          error stop
        End select
        i = 0
        do nc1 = 1, N_Fam
          i = i + L_Fam(nc1)
        enddo
        if( nc /= i ) then
          Write(error_unit,*) "Error in Predefined_checkerboard", nc, N_Fam, L_Fam
          error stop
        endif
           
      End Subroutine Predefined_checkerboard
      
      Subroutine Checkeboard_Add_subfam(checker_in, L_Fam0, N_unit, List0, N_orb, N_sub, &
                                        List2, Invlist2, checker_inter, L_Fam_inter, checker_intra)

        Implicit none
        
        Integer, Intent(in)               :: checker_in(:,:,:), L_Fam0(:), List0(:,:)
        Integer, Intent(in)               :: N_unit, N_orb, N_sub
        Integer, Intent(out), allocatable :: List2(:,:), Invlist2(:,:,:)
        Integer, Intent(out), allocatable :: checker_inter(:,:,:)
        Integer, Intent(out), allocatable :: L_Fam_inter(:)
        Integer, Intent(out), allocatable :: checker_intra(:,:,:)

        !Local
        Integer :: N_Fam0, N_Fam, L_Fam_max, B_intra, L_intra
        Integer :: n, n2, l, l2, d, il1, il2, nc1, nc2, I1, I2, no1, no2, nb, nc, I, no
        
        
        call Latt_add_orbitals(N_unit, N_orb, N_sub, List2, Invlist2)
        
        
        !TODO: Optimize checker_intra for N_sub>4
        select case(N_sub)
        case(4)
          B_intra = 3
          L_intra = N_unit*N_orb*2
          allocate( checker_intra(2, L_intra, B_intra) )
          nc = 0
          do I = 1, N_unit
            do no = 1,N_orb
              nc = nc+1
              nb = 1
              checker_intra(1, nc, nb) = Invlist2(I,no,1)
              checker_intra(2, nc, nb) = Invlist2(I,no,2)
              nb = 2
              checker_intra(1, nc, nb) = Invlist2(I,no,1)
              checker_intra(2, nc, nb) = Invlist2(I,no,3)
              nb = 3
              checker_intra(1, nc, nb) = Invlist2(I,no,1)
              checker_intra(2, nc, nb) = Invlist2(I,no,4)
              nc = nc+1
              nb = 1
              checker_intra(1, nc, nb) = Invlist2(I,no,3)
              checker_intra(2, nc, nb) = Invlist2(I,no,4)
              nb = 2
              checker_intra(1, nc, nb) = Invlist2(I,no,2)
              checker_intra(2, nc, nb) = Invlist2(I,no,4)
              nb = 3
              checker_intra(1, nc, nb) = Invlist2(I,no,2)
              checker_intra(2, nc, nb) = Invlist2(I,no,3)
            enddo
          enddo
              
        case default
          B_intra = (N_sub*(N_sub-1))/2
          allocate( checker_intra(2, N_unit*N_orb, B_intra) )
          nb = 0
          do n = 1,N_sub-1
            do d = 1, N_sub-n
              nb = nb+1
              nc = 0
              do I = 1, N_unit
                do no = 1,N_orb
                  nc = nc+1
                  nc1 = Invlist2(I,no,n)
                  nc2 = Invlist2(I,no,n+d)
                  checker_intra(1, nc, nb) = nc1
                  checker_intra(2, nc, nb) = nc2
                enddo
              enddo
            enddo
          enddo
          if ( nb /= B_intra ) then
            write(error_unit,*) "Error: nb /= B_intra", nb, B_intra
            error stop
          endif
        end select
        
        
        !L_Fam0 = size(checker_in, 2)
        N_Fam0 = size(checker_in, 3)
        L_Fam_max  = maxval(L_Fam0) * N_sub
        N_Fam  = N_Fam0 * N_sub
        allocate( L_Fam_inter(N_Fam), checker_inter(2, L_Fam_max, N_fam) )
        
        n2 = 0
        do n = 1,N_Fam0
          do d = 0, N_sub-1
            n2 = n2+1
            l2 = 0
            do l = 1,L_Fam0(n)
              do il1 = 1, N_sub
                l2 = l2+1
                nc1 = checker_in(1, l, n)
                nc2 = checker_in(2, l, n)
                il2 = il1+d
                if (il2 > N_sub) il2 = il2-N_sub
                
                I1  = List0(nc1,1)
                no1 = List0(nc1,2)
                I2  = List0(nc2,1)
                no2 = List0(nc2,2)
                
                checker_inter(1, l2, n2) = Invlist2(I1,no1,il1)
                checker_inter(2, l2, n2) = Invlist2(I2,no2,il2)
              enddo
            enddo
            L_Fam_inter(n2) = l2
          enddo
        enddo
        
      End Subroutine Checkeboard_Add_subfam
      
    end Module Predefined_structures
