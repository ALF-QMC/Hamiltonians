!  Copyright (C) 2022 The ALF project
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
!       http://alf.physik.uni-wuerzburg.de
!
!     - We require the preservation of the above copyright notice and this license in all original files.
!
!     - We prohibit the misrepresentation of the origin of the original source files. To obtain
!       the original source files please visit the homepage http://alf.physik.uni-wuerzburg.de .
!
!     - If you make substantial changes to the program we require you to either consider contributing
!       to the ALF project or to mark your material in a reasonable way as different from the original version


!--------------------------------------------------------------------
!> @author
!> ALF-project
!>
!> @brief
!> This module defines the  Hamiltonian and observables.  Here, we have included a 
!> Quantum-Spin-Hall Model with SU(N) color index.
!>
!> @details
!> The public variables of this module are the following
!>
!> @param [public] OP_V
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> List of operators of type=1,2 and 3 describing the sequence of interactions on a time slice.
!> The first index runs over this sequence. The second corresponds to the flavor index.  \endverbatim
!>
!> @param [public] OP_T
!> \verbatim
!> Type (Operator), dimension(:,:), allocatable
!> Sequence of  operators  accounting for the  hopping on a  time slice. This can include  various
!> checkerboard decompositions. The first index runs over this sequence. The second corresponds to
!> the flavor index. \endverbatim
!> *  The progagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n}  \f$.  That is
!> first the hopping and then the potential energy.
!>
!>@param [public] WF_L
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Left trial wave function.  \endverbatim
!>
!> @param [public] WF_R
!> \verbatim Type (WaveFunction), dimension(:),   allocatable
!> Right trial wave function.   For both wave functions the index runs over the flavor index. \endverbatim
!>
!> @param [public]  nsigma
!> \verbatim Type(Fields)
!> Contains all auxiliary fields in the variable f(:,:). The first index runs through the operator
!> sequence. The second through the time slices.   \endverbatim
!
!> @param [public]  Ndim
!> \verbatim Integer
!> Total number of orbitals. e.g. # unit cells * # orbitals per unit cell.  \endverbatim
!
!> @param [public]  N_FL
!> \verbatim Integer
!> # of flavors.  Propagation is block diagonal in flavors.  \endverbatim
!
!> @param [public]  N_SUN
!> \verbatim Integer
!> # of colors.  Propagation is color independent.  \endverbatim
!>
!> @param [public] Ltrot
!> \verbatim Integer
!> Available measurment interval in units of Delta Tau. \endverbatim
!>
!> @param [public] Thtrot
!>  \verbatim Integer
!> Effective projection parameter in units of Delta Tau.  (Only relevant if projective option is turned on) \endverbatim
!>
!> @param [public] Projector
!> \verbatim Logical
!> Flag for projector. If true then the total number of time slices will correspond to Ltrot + 2*Thtrot \endverbatim
!>
!> @param [public] Group_Comm
!> \verbatim Integer
!> Defines MPI communicator  \endverbatim
!
!> @param [public] Symm
!> \verbatim Logical  \endverbatim
!> If set to true then the green functions will be symmetrized
!> before being  sent to the Obser, ObserT subroutines.
!> In particular, the transformation,  \f$ \tilde{G} =  e^{-\Delta \tau T /2 } G e^{\Delta \tau T /2 } \f$
!> will be carried out  and \f$ \tilde{G} \f$  will be sent to the Obser and ObserT subroutines.  Note that
!> if you want to use this  feature, then you have to be sure the hopping and interaction terms are decomposed
!> symmetrically. If Symm is true, the propagation reads:
!> \f$ \prod_{\tau} \; \;  \prod_{n=N_T}^{1}e^{T_n/2} \prod_{n=1}^{N_V}e^{V_n(\tau)}  \prod_{n=1}^{N_T}e^{T_n/2}  \f$
!>
!>
!> You still have to add some docu for the other private variables in this module.
!>
!--------------------------------------------------------------------

    submodule (Hamiltonian_main) ham_QSH_smod

    Use Operator_mod
    Use WaveFunction_mod
    Use Lattices_v3
    Use MyMats
    Use Random_Wrap
    Use Files_mod
    Use Matrix
    Use Observables
    Use Fields_mod
    Use Predefined_Hoppings
    Use LRC_Mod

    Implicit none
    
    type, extends(ham_base) :: ham_QSH
    contains
      ! Set Hamiltonian-specific procedures
      procedure, nopass :: Ham_Set
      procedure, nopass :: Alloc_obs
      procedure, nopass :: Obser
      procedure, nopass :: ObserT

#ifdef HDF5
      procedure, nopass :: write_parameters_hdf5
#endif
    end type ham_QSH

    !#PARAMETERS START# VAR_lattice
    Character (len=64) :: Model = ''  ! Value irrelevant
    Character (len=64) :: Lattice_type = 'Bilayer_honeycomb'  ! Possible Values: 'Bilayer_honeycomb'
    Integer            :: L1 = 6   ! Length in direction a_1
    Integer            :: L2 = 6   ! Length in direction a_2
    !#PARAMETERS END#

    !#PARAMETERS START# VAR_Model_Generic
    !Integer              :: N_SUN        = 2        ! Number of colors, globally defined in Hamiltonian_main_mod.F90
    !Integer              :: N_FL         = 1        ! Number of flavors
    real(Kind=Kind(0.d0)) :: Phi_X        = 0.d0     ! Twist along the L_1 direction, in units of the flux quanta
    real(Kind=Kind(0.d0)) :: Phi_Y        = 0.d0     ! Twist along the L_2 direction, in units of the flux quanta
    logical               :: Bulk         = .true.   ! Twist as a vector potential (.T.), or at the boundary (.F.)
    Integer               :: N_Phi        = 0        ! Total number of flux quanta traversing the lattice
    real(Kind=Kind(0.d0)) :: Dtau         = 0.1d0    ! Thereby Ltrot=Beta/dtau
    real(Kind=Kind(0.d0)) :: Beta         = 12.d0     ! Inverse temperature
    logical               :: Checkerboard = .true.   ! Whether checkerboard decomposition is used
    !logical              :: Symm         = .true.   ! Whether symmetrization takes place
    !logical              :: Projector    = .false.  ! Whether the projective algorithm is used
    real(Kind=Kind(0.d0)) :: Theta        = 10.d0    ! Projection parameter
    !#PARAMETERS END#

    !#PARAMETERS START# VAR_QSH
    real(Kind=Kind(0.d0)) :: Ham_T      = 1.d0      ! Hopping parameter
    real(Kind=Kind(0.d0)) :: Ham_T2     = 1.d0      ! Hopping parmaeter layer 2
    real(Kind=Kind(0.d0)) :: Ham_chem   = 0.d0      ! Chemical potential
    real(Kind=Kind(0.d0)) :: Ham_lambda = 0.01d0    ! Interaction 
    !#PARAMETERS END#

    Type (Lattice),   target :: Latt
    Type (Unit_cell), target :: Latt_unit        ! Regular unit cell of bilayer honeycomb
    Type (Unit_cell), target :: Latt_Unit_Hex    ! Full hexagon unit cell
    Type (Unit_cell), target :: Latt_Unit_Obs    ! Unit cell of single layer honeycomb for observables
    Type (Unit_cell), target :: Latt_Unit_nnBond ! Three nearest neighbor bonds
    Type (Unit_cell), target :: Latt_Unit_ABC    ! Unit cell for 3 groups A, B and C

    Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
    Integer, allocatable :: List(:,:), Invlist(:,:)                     ! For orbital structure of Unit cell
    Integer, allocatable :: List_Obs(:,:), Invlist_Obs(:,:)             ! For "physical" unit cell of honeycomb lattice
    Integer, allocatable :: List_Bonds(:,:), Invlist_Bonds(:,:,:)       ! For next-nearest neighbor bonds inside hexagon
    Integer, allocatable :: List_nnBonds(:,:), Invlist_nnBonds(:,:,:)   ! For nearest neighbor bonds
    Integer, allocatable :: List_A(:), List_B(:), List_C(:)             ! For different hexagon groups
    Integer, allocatable :: Hexagon_List(:,:)                           ! List for filling Op_V%P containing the hexagon orbitals

    real (Kind=Kind(0.d0)) :: Ham_Tperp

  contains
    
    module Subroutine Ham_Alloc_QSH
      allocate(ham_QSH::ham)
    end Subroutine Ham_Alloc_QSH

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_QSH_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
    Subroutine Ham_Set

#if defined (MPI) || defined(TEMPERING)
        Use mpi
#endif
        Implicit none

        integer                :: ierr, nf, unit_info
        Character (len=64)     :: file_info


#ifdef MPI
        Integer        :: Isize, Irank, irank_g, isize_g, igroup
        Integer        :: STATUS(MPI_STATUS_SIZE)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        call MPI_Comm_rank(Group_Comm, irank_g, ierr)
        call MPI_Comm_size(Group_Comm, isize_g, ierr)
        igroup           = irank/isize_g
#endif

        ! From dynamically generated file "Hamiltonian_QSH_read_write_parameters.F90"
        call read_parameters()

        Ltrot = nint(beta/dtau)
        if (Projector) Thtrot = nint(theta/dtau)
        Ltrot = Ltrot+2*Thtrot

        Ham_Tperp = 0.d0

        ! Setup the Bravais lattice
        call Ham_Latt

        ! Setup the hopping / single-particle part
        call Ham_Hop

        ! Setup the interaction.
        call Ham_V

        ! Setup the trial wave function, in case of a projector approach
        if (Projector) Call Ham_Trial(file_info)
#ifdef MPI
        If (Irank_g == 0) then
#endif
           File_info = "info"
#if defined(TEMPERING)
           write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
           Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      :  QSH'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# unit cells  : ', Latt%N 
             Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
             Write(unit_info,*) 'Checkerboard  : ', Checkerboard
             Write(unit_info,*) 'Symm. decomp  : ', Symm
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
             else
                Write(unit_info,*) 'Finite temperture version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ',   N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 't             : ', Ham_T
             Write(unit_info,*) 'Ham_lambda    : ', Ham_lambda
             Write(unit_info,*) 'Ham_chem      : ', Ham_chem

             if (Projector) then
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             endif
           Close(unit_info)
#ifdef MPI
        Endif
#endif
      end Subroutine Ham_Set


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the  Lattice
!--------------------------------------------------------------------
      Subroutine Ham_Latt

        Use Predefined_Lattices

        Implicit none

        Integer :: nc, no, i, j, i1, i2, I_nn1, I_nn2, I_nn3, Is_up, Is_dn, Js_up, Js_dn
        ! Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)

        ! Use predefined stuctures or set your own lattice.
        If ( Lattice_type .ne. "Bilayer_honeycomb" )  then
           write(6,*)  ' Program only runs for bilayer honeycomb '
           stop
        elseif ((mod(L1,3) .ne. 0) .or. (mod(L2,3) .ne. 0)) then
           write(6,*)  ' L1 and L2 must be multiples of 3 '
           stop
        else
           call Predefined_Latt(Lattice_type, L1,L2,Ndim, List,Invlist,Latt,Latt_Unit)
        endif

        ! print * , 'i   ', 'i_1=Latt%List(i,1)  ', 'i_2=Latt%List(i,2)  ', 'Invlist(i_1,i_2)'
        ! do i = 1,Latt%N
        !    print * , i, Latt%List(i,1), Latt%List(i,2), Latt%invlist(Latt%List(i,1),Latt%List(i,2))
        ! enddo

        ! Fill lists for unit cells of groups A, B and C 
        Allocate (List_A(Latt%N/3), List_B(Latt%N/3), List_C(Latt%N/3))
        i1 = 1
        do i = 1,L2 
           i2 = i1
           do j = 1,L1/3               
              List_A(i+L1*(j-1)) = i2
              List_B(i+L1*(j-1)) = Latt%nnlist(i2,1,0)
              List_C(i+L1*(j-1)) = Latt%nnlist(i2,0,1)
              i2 = Latt%nnlist(Latt%nnlist(Latt%nnlist(i2,1,0),1,0),1,0)
           enddo
           i1 = Latt%nnlist(i1,1,1)
        enddo
        ! print *, List_A
        ! print *, List_B 
        ! print *, List_C


        ! Convention:
        !       o 3
        !      /|\
        !     / | \
        !    / / \ \
        ! 4 o v   ^ o 6
        !   |/+   +\|
        ! 1 o--->---o 2
        !    \  +  /
        !     \   /
        !      \ /
        !       o
        !       5
        !
        !       o 3
        !      / \
        !     /   \
        !    /  +  \
        ! 4 o---<---o 6
        !   |\+   +/|
        ! 1 o v   ^ o 2
        !    \ \ / /
        !     \ | /
        !      \|/
        !       o
        !       5
        Allocate (Hexagon_List(Latt%N,12))
        do i = 1,Latt%N 
           Hexagon_List(i,1)  = Invlist(I,1)    !  Orbital 1 spin up --> Site 1 of hexagon
           Hexagon_List(i,2)  = Invlist(I,3)    !  Orbital 1 spin down
           Hexagon_List(i,3)  = Invlist(Latt%nnlist(I,1,0),1)  !  Site 2 spin up
           Hexagon_List(i,4)  = Invlist(Latt%nnlist(I,1,0),3)
           Hexagon_List(i,5)  = Invlist(Latt%nnlist(I,0,1),1)  !  Site 3 spin up
           Hexagon_List(i,6)  = Invlist(Latt%nnlist(I,0,1),3)
           Hexagon_List(i,7)  = Invlist(I,2)    !  Orbital 2 spin up --> Site 4 of hexagon
           Hexagon_List(i,8)  = Invlist(I,4)    !  Orbital 2 spin down
           Hexagon_List(i,9)  = Invlist(Latt%nnlist(I,1,-1),2) !  Site 5 spin up
           Hexagon_List(i,10) = Invlist(Latt%nnlist(I,1,-1),4)
           Hexagon_List(i,11) = Invlist(Latt%nnlist(I,1,0),2)  !  Site 6 spin up
           Hexagon_List(i,12) = Invlist(Latt%nnlist(I,1,0),4)
        enddo

        Allocate (List_Bonds(Latt%N*6,2), Invlist_Bonds(Latt%N,6,4))
        nc = 0
        Do I = 1,Latt%N
           Do no = 1,6
              nc = nc + 1

              List_Bonds(nc,1) = I
              List_Bonds(nc,2) = no

              I_nn1 = Latt%nnlist(I,1,0); I_nn2 = Latt%nnlist(I,0,1)
              I_nn3 = Latt%nnlist(I,1,-1)

              select case (no)
              case (1)   ! Bond 1 -> 2 (Is = 1, Js = 2)
                  Is_up = Invlist(I    ,1); Is_dn = Invlist(I    ,3)
                  Js_up = Invlist(I_nn1,1); Js_dn = Invlist(I_nn1,3)
              case (2)   ! Bond 2 -> 3
                  Is_up = Invlist(I_nn1,1); Is_dn = Invlist(I_nn1,3)
                  Js_up = Invlist(I_nn2,1); Js_dn = Invlist(I_nn2,3)
              case (3)   ! Bond 3 -> 1
                  Is_up = Invlist(I_nn2,1); Is_dn = Invlist(I_nn2,3)
                  Js_up = Invlist(I    ,1); Js_dn = Invlist(I    ,3)
              case (4)   ! Bond 4 -> 5
                  Is_up = Invlist(I    ,2); Is_dn = Invlist(I    ,4)
                  Js_up = Invlist(I_nn3,2); Js_dn = Invlist(I_nn3,4)
              case (5)   ! Bond 5 -> 6
                  Is_up = Invlist(I_nn3,2); Is_dn = Invlist(I_nn3,4)
                  Js_up = Invlist(I_nn1,2); Js_dn = Invlist(I_nn1,4)
              case (6)   ! Bond 6 -> 4
                  Is_up = Invlist(I_nn1,2); Is_dn = Invlist(I_nn1,4)
                  Js_up = Invlist(I    ,2); Js_dn = Invlist(I    ,4)
              case default
                 Write(6,*) ' Error in orbital setting '
              end select
              Invlist_Bonds(I,no,1) =  Is_up
              Invlist_Bonds(I,no,2) =  Is_dn
              Invlist_Bonds(I,no,3) =  Js_up
              Invlist_Bonds(I,no,4) =  Js_dn
           Enddo
        Enddo

        Allocate (List_nnBonds(Latt%N*3,2), Invlist_nnBonds(Latt%N,3,4))
        nc = 0
        Do I = 1,Latt%N
           Do no = 1,3
              nc = nc + 1

              List_nnBonds(nc,1) = I
              List_nnBonds(nc,2) = no

              I_nn1 = Latt%nnlist(I,1,-1)
              I_nn2 = Latt%nnlist(I,0,-1)

              select case (no)
              case (1)   ! Bond 1 -> 4 (Is = 1, Js = 4)
                  Is_up = Invlist(I    ,1); Is_dn = Invlist(I    ,3)
                  Js_up = Invlist(I    ,2); Js_dn = Invlist(I    ,4)
              case (2)   ! Bond 1 -> 5
                  Is_up = Invlist(I    ,1); Is_dn = Invlist(I    ,3)
                  Js_up = Invlist(I_nn1,2); Js_dn = Invlist(I_nn1,4)
              case (3)   ! Bond 3 -> 4
                  Is_up = Invlist(I    ,1); Is_dn = Invlist(I    ,3)
                  Js_up = Invlist(I_nn2,2); Js_dn = Invlist(I_nn2,4)
              case default
                 Write(6,*) ' Error in orbital setting '
              end select
              Invlist_nnBonds(I,no,1) =  Is_up
              Invlist_nnBonds(I,no,2) =  Is_dn
              Invlist_nnBonds(I,no,3) =  Js_up
              Invlist_nnBonds(I,no,4) =  Js_dn
           Enddo
        Enddo

        ! Latt Unit describing position of next-nearest neighbor bonds (i.e. involved sites)
        Latt_Unit_Hex%Norb    = 6
        Latt_Unit_Hex%N_coord = 6 ! Each site has 6 next-nearest neighbors
        Allocate (Latt_Unit_Hex%Orb_pos_p(Latt_Unit_Hex%Norb,2))
        do nc = 1,2
           Latt_Unit_Hex%Orb_pos_p(1,nc) = 0.d0
           Latt_Unit_Hex%Orb_pos_p(2,nc) = Latt%a1_p(nc)
           Latt_Unit_Hex%Orb_pos_p(3,nc) = Latt%a2_p(nc)
           Latt_Unit_Hex%Orb_pos_p(4,nc) = (Latt%a2_p(nc) - 0.5D0*Latt%a1_p(nc) ) * 2.D0/3.D0
           Latt_Unit_Hex%Orb_pos_p(5,nc) = (2.D0*Latt%a1_p(nc) - Latt%a2_p(nc) ) * 1.D0/3.D0
           Latt_Unit_Hex%Orb_pos_p(6,nc) = (Latt%a2_p(nc) + Latt%a1_p(nc) ) * 2.D0/3.D0
        enddo


        ! Latt_Unit describing position of nearest neighbor bonds
        Latt_Unit_nnBond%Norb    = 3
        Latt_Unit_nnBond%N_coord = 3 ! each site has 3 nearest neighbors
        Allocate (Latt_Unit_nnBond%Orb_pos_p(3,3))
        do nc = 1,2
           Latt_Unit_nnBond%Orb_pos_p(1,nc) = (Latt%a2_p(nc) - 0.5D0*Latt%a1_p(nc) ) * 1.D0/3.D0
           Latt_Unit_nnBond%Orb_pos_p(2,nc) = (2.D0*Latt%a1_p(nc) - Latt%a2_p(nc) ) * 1.D0/6.D0
           Latt_Unit_nnBond%Orb_pos_p(3,nc) = -(Latt%a2_p(nc) + Latt%a1_p(nc) ) * 1.D0/6.D0
        enddo

        Latt_Unit_Obs%Norb    = 2
        Latt_Unit_Obs%N_coord = 3
        Allocate (Latt_Unit_Obs%Orb_pos_p(2,2))
        Latt_Unit_Obs%Orb_pos_p(1,:) = 0.d0
        Latt_Unit_Obs%Orb_pos_p(2,:) = (Latt%a2_p(:) - 0.5D0*Latt%a1_p(:) ) * 2.D0/3.D0

        Allocate (List_Obs(Latt%N*Latt_Unit_Obs%Norb,2), Invlist_Obs(Latt%N,Latt_Unit_Obs%Norb))
         nc = 0
         Do I = 1,Latt%N
            Do no = 1,Latt_Unit_Obs%Norb
               nc = nc + 1
               List_Obs(nc,1) = I
               List_Obs(nc,2) = no 
               Invlist_Obs(I,no) = nc
            Enddo
         Enddo

       end Subroutine Ham_Latt


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets  the Hopping
!--------------------------------------------------------------------
       Subroutine Ham_Hop

        Implicit none

        Real (Kind=Kind(0.d0) ) ::  Ham_Lambda = 0.d0

        Real (Kind=Kind(0.d0) ), allocatable :: Ham_T_vec(:), Ham_Tperp_vec(:), Ham_Chem_vec(:), Phi_X_vec(:), Phi_Y_vec(:),&
           &                                  Ham_T2_vec(:),  Ham_Lambda_vec(:)
        Integer, allocatable ::   N_Phi_vec(:)

        ! Use predefined stuctures or set your own hopping
        Integer :: n,nth

        Allocate (Ham_T_vec(N_FL), Ham_T2_vec(N_FL), Ham_Tperp_vec(N_FL), Ham_Chem_vec(N_FL), Phi_X_vec(N_FL), Phi_Y_vec(N_FL),&
           &                                   N_Phi_vec(N_FL), Ham_Lambda_vec(N_FL) )

        ! Here we consider no N_FL  dependence of the hopping parameters.
        Ham_T_vec      = Ham_T
        Ham_Tperp_vec  = Ham_Tperp
        Ham_Chem_vec   = Ham_Chem
        Phi_X_vec      = Phi_X
        Phi_Y_vec      = Phi_Y
        Ham_T2_vec     = Ham_T2
        Ham_Lambda_vec = Ham_Lambda
        N_Phi_vec      = N_Phi

        Select case (Lattice_type)
        Case ("Square")
           Call  Set_Default_hopping_parameters_square(Hopping_Matrix,Ham_T_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
              &                                      Bulk, N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
        Case ("N_leg_ladder")
           Call  Set_Default_hopping_parameters_n_leg_ladder(Hopping_Matrix, Ham_T_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, &
              &                                            Phi_Y_vec, Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
        Case ("Honeycomb")
           Ham_Lambda = 0.d0
           Call  Set_Default_hopping_parameters_honeycomb(Hopping_Matrix, Ham_T_vec, Ham_Lambda_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
              &                                         Bulk,  N_Phi_vec, N_FL, List, Invlist, Latt, Latt_unit )
        Case ("Bilayer_square")
           Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
              &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
              &                                              List, Invlist, Latt, Latt_unit )

        Case ("Bilayer_honeycomb")
           Call  Set_Default_hopping_parameters_Bilayer_honeycomb(Hopping_Matrix,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
              &                                                 Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
              &                                                 List, Invlist, Latt, Latt_unit )

        end Select

        Call  Predefined_Hoppings_set_OPT(Hopping_Matrix,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_T )

        Deallocate (Ham_T_vec, Ham_T2_vec, Ham_Tperp_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec, &
           &                                   N_Phi_vec,  Ham_Lambda_vec )


    end Subroutine Ham_Hop
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the trial wave function
!--------------------------------------------------------------------
    Subroutine Ham_Trial(file_info)


#if defined (MPI) || defined(TEMPERING)
      Use mpi
#endif
      Use Predefined_Trial

      Implicit none
      Character (len=64), intent(in)  :: file_info


      Integer :: N_part, nf
#ifdef MPI
      Integer        :: Isize, Irank, irank_g, isize_g, igroup, ierr
      Integer        :: STATUS(MPI_STATUS_SIZE)
      Type(Operator),  dimension(:,:), allocatable  :: OP_tmp
      Type (Lattice)                                :: Latt_Kekule
      Real (Kind=Kind(0.d0))  :: A1_p(2), A2_p(2), L1_p(2), L2_p(2), x_p(2),x1_p(2), hop(3)
      Real (Kind=Kind(0.d0))  :: delta = 0.01
      Integer :: N, I, I1, I2, nc, nc1, IK_u, I_u, J1, ns
      Complex (Kind=Kind(0.d0)) :: Z_norm

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
      call MPI_Comm_rank(Group_Comm, irank_g, ierr)
      call MPI_Comm_size(Group_Comm, isize_g, ierr)
      igroup           = irank/isize_g
#endif
      ! Use predefined stuctures or set your own Trial  wave function
      N_part = Ndim/2
      ! Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
      !    &                            N_part, N_FL,  WF_L, WF_R)

      ! Ham_T_vec     = 1.d0
      ! Ham_T2_vec    = 1.d0
      ! Ham_Tperp_vec = 0.d0
      ! Phi_X_vec     = 0.00
      ! Call  Set_Default_hopping_parameters_Bilayer_square(Hopping_Matrix_tmp,Ham_T_vec,Ham_T2_vec,Ham_Tperp_vec, Ham_Chem_vec, &
      !    &                                              Phi_X_vec, Phi_Y_vec, Bulk,  N_Phi_vec, N_FL,&
      !    &                                              List, Invlist, Latt, Latt_unit )
      ! Call  Predefined_Hoppings_set_OPT(Hopping_Matrix_tmp,List,Invlist,Latt,  Latt_unit,  Dtau, Checkerboard, Symm, OP_tmp )
      !  Kekule Mass term to avoid  degeneracy at half-filling.

      Allocate(Op_Tmp(1,N_FL))
      do n = 1,N_FL
         Call Op_make(Op_Tmp(1,n),Ndim)
      enddo
      Allocate(WF_L(N_FL),WF_R(N_FL))
        do n=1,N_FL
           Call WF_alloc(WF_L(n),Ndim,N_part)
           Call WF_alloc(WF_R(n),Ndim,N_part)
        enddo

      A1_p = 2.d0 * Latt%a1_p  - Latt%a2_p
      A2_p =        Latt%a1_p  + Latt%a2_p
      L1_p = Latt%L1_p
      L2_p = Latt%L2_p
      Call Make_Lattice( L1_p, L2_p, A1_p,  A2_p, Latt_Kekule)

      DO I = 1, Latt_Kekule%N
         x_p = dble(Latt_Kekule%list(I,1))*Latt_Kekule%a1_p + dble(Latt_Kekule%list(I,2))*Latt_Kekule%a2_p
         IK_u   = Inv_R(x_p,Latt)
         do nc  = 1, 3
            select case (nc)
            case (1)
               I_u    =  IK_u
               hop(1) =  1.d0 + delta
               hop(2) =  1.d0 - delta
               hop(3) =  1.d0
            case (2)
               I_u    = Latt%nnlist(IK_u,0,1)
               hop(1) =  1.d0
               hop(2) =  1.d0 + delta
               hop(3) =  1.d0 - delta
            case (3)
               I_u     = Latt%nnlist(IK_u,1,0)
               hop(1) =  1.d0 - delta
               hop(2) =  1.d0
               hop(3) =  1.d0 + delta
            end select
            x_p = dble(Latt%list(I_u,1))*Latt%a1_p + dble(Latt%list(I_u,2))*Latt%a2_p
            do ns = 0,1
               I1 = invlist(I_u,1+2*ns)
               do nc1 = 1,3
                  select case (nc1)
                  case (1)
                     J1 = invlist(I_u,2+2*ns)
                  case (2)
                     J1 = invlist(Latt%nnlist(I_u,1,-1),2+2*ns)
                  case (3)
                     J1 = invlist(Latt%nnlist(I_u,0,-1),2+2*ns)
                  end select

                  do n = 1,N_FL
                     Op_Tmp(1,n)%O(I1,J1) =   cmplx( - hop(nc1),    0.d0, kind(0.D0))
                     Op_Tmp(1,n)%O(J1,I1) =   cmplx( - hop(nc1),    0.d0, kind(0.D0))
                  enddo
               enddo
            enddo
         enddo
      Enddo
      do n = 1,N_FL
         Do I = 1,Ndim
            Op_Tmp(1,n)%P(i) = i
         Enddo
         Op_Tmp(1,n)%g    = cmplx(1.d0, 0.d0,kind(0.d0))
         Op_Tmp(1,n)%alpha= cmplx(0.d0,0.d0, kind(0.D0))
         Call Op_set(Op_Tmp(1,n))
      Enddo
      
      Do nf = 1,N_FL
         Call Diag(Op_tmp(1,nf)%O,Op_tmp(1,nf)%U,Op_tmp(1,nf)%E)
         do I2=1,N_part
            do I1=1,Ndim
               WF_L(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
               WF_R(nf)%P(I1,I2)=Op_tmp(1,nf)%U(I1,I2)
            enddo
         enddo
         WF_L(nf)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
         WF_R(nf)%Degen = Op_tmp(1,nf)%E(N_part+1) - Op_tmp(1,nf)%E(N_part)
      enddo

      Do nf = 1,N_FL
         Call WF_overlap(WF_L(nf), WF_R(nf), Z_norm)
         !Write(6,*) " Z_norm ", Z_norm
      enddo

      Do nf = 1,N_FL
         Call Op_clear(OP_tmp(1,nf),Ndim)
      enddo
      Deallocate (OP_tmp)
      !Call Predefined_hoppings_clear(Hopping_Matrix_tmp)
      !Deallocate (Ham_T_vec, Ham_Tperp_vec, Ham_T2_vec, Ham_Chem_vec, Phi_X_vec, Phi_Y_vec,  N_Phi_vec )

#ifdef MPI
      If (Irank_g == 0) then
#endif
         OPEN(Unit = 50,file=file_info,status="unknown",position="append")
         Do nf = 1,N_FL
            Write(50,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
            Write(50,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
         enddo
         close(50)
#ifdef MPI
      endif
#endif

    end Subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
    Subroutine Ham_V

     Use Predefined_Int
     Implicit none

     Integer :: I, is, N_op, nc, no, nf, I1up, I1dn, I2up, I2dn, I3up, I3dn, I4up, I4dn, I5up, I5dn, I6up, I6dn
     Integer :: nb, ni, n1, n2, n3, n4
     Complex(Kind=Kind(0.d0))          :: Z1, Z2
     Real(Kind=Kind(0.d0))             :: Zero =  1.0E-6

     
     ! Convention:
     !       o 3
     !      /|\
     !     / | \
     !    / / \ \
     ! 4 o ^   ^ o 6
     !   |/-   +\|
     ! 1 o--->---o 2
     !    \  +  /
     !     \   /
     !      \ /
     !       o
     !       5
     !             / 0  1 -1 \
     ! \nu_{ij} = | -1  0  1  | = \epsilon_{ij} for i,j={1,2,3},{4,5,6} respectively
     !             \ 1 -1  0 /
     !
     !
     !                       /      0      i\sigma^a   -i\sigma^a       0           0           0     \
     !                      / -i\sigma^a        0       i\sigma^a       0           0           0      \
     ! --> Op_V%O(i,j) =   |   i\sigma^a  -i\sigma^a        0           0           0           0       | for a={x,y,z}
     !                     |        0           0           0           0       i\sigma^a   -i\sigma^a  | and i,j={1,...,12}
     !                      \       0           0           0      -i\sigma^a       0        i\sigma^a /
     !                       \      0           0           0       i\sigma^a  -i\sigma^a       0     /
     !

     N_op = 0
     if (abs(Ham_lambda) > Zero ) Then
        ! Number of interaction terms:
        ! H_\lambda decomposes into 3*3 (sigma,hexagon groups) parts with Latt%N/3 terms
        ! -> Total of (9*2-1)*Latt%N/3 terms 
        N_op = 17 * Latt%N/3 
     endif
     Allocate(Op_V(N_op,N_FL))

     if ( abs(Ham_lambda)  > Zero ) then
        do nf = 1,N_FL
           do i  = 1, N_op
              Call Op_make(Op_V(i,nf), 12)
           enddo
        enddo

        Do nf = 1,N_FL
           nc = 0
           
           !  Hexagon group A
           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)
              
              ! Write Op_V(nc,nf)%P(i)
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo

              ! Write Op_V(nc,nf)%O(i,j)
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo
              

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_z  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+3 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2+1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = cmplx( 0.d0,-dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  Hexagon group B
           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)
              
              ! Write Op_V(nc,nf)%P(i)
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
              
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_z  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+3 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2+1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = cmplx( 0.d0,-dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  Hexagon group C
           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_C(is)
              
              ! Write Op_V(nc,nf)%P(i)
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
              
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_C(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_z  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_C(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+3 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2+1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = cmplx( 0.d0,-dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN)), 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_C(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_C(is)
              
              ! Write Op_V(nc,nf)%P(i)
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
              
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  Hexagon group B
           !  sigma_z  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+3 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2+1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = cmplx( 0.d0,-dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5D0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_B(is)
              
              ! Write Op_V(nc,nf)%P(i)
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
              
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo
           !  Hexagon group A
           !  sigma_z  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+3 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2+1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = cmplx( 0.d0,-dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5D0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_y  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)

              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
                 
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Z2 = cmplx(-dble(1-(ni/2)*2+(ni/3)*2), 0.d0, kind(0.D0))
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

           !  sigma_x  interaction
           do is = 1,Latt%N/3
              nc = nc + 1
              I = List_A(is)
              
              do no = 1,12
                 Op_V(nc,nf)%P(no) = Hexagon_List(I,no)
              enddo
              
              do nb = 0,1
                 do ni = 1,3
                    n1 = (ni/3)*2+1 + nb*6     ! row of first element
                    n2 = ((ni+1)/3)*2+4 + nb*6 ! column of first element
                    n3 = (n1+1)                ! row of second element
                    n4 = (n2-1)                ! column of second element
                    Z1 = cmplx( 0.d0, dble(1-(ni/2)*2+(ni/3)*2), kind(0.D0))
                    Z2 = Z1
                    Op_V(nc,nf)%O(n1,n2)=Z1 ; Op_V(nc,nf)%O(n3,n4)=Z2
                    Op_V(nc,nf)%O(n2,n1)=Conjg(Z1) ; Op_V(nc,nf)%O(n4,n3)=Conjg(Z2)
                 enddo
              enddo

              Op_V(nc,nf)%g      = SQRT(CMPLX(DTAU*Ham_lambda/(DBLE(N_SUN))*0.5d0, 0.D0, kind(0.D0)))
              Op_V(nc,nf)%alpha  = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,nf)%type   = 2
              Call Op_set( Op_V(nc,nf) )
           enddo

        Enddo
     endif

      
    end Subroutine Ham_V




!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
    Subroutine  Alloc_obs(Ltau)

        Implicit none
        !>  Ltau=1 if time displaced correlations are considered.
        Integer, Intent(In) :: Ltau
        Integer    ::  i, N, Nt
        Character (len=64) ::  Filename
        Character (len=2)  ::  Channel


        ! Scalar observables
        Allocate ( Obs_scal(5) )
        Do I = 1,Size(Obs_scal,1)
           select case (I)
           case (1)
              N = 1;   Filename ="Kin"
           case (2)
              N = 1;   Filename ="Pot"
           case (3)
              N = 1;   Filename ="Part"
           case (4)
              N = 1;   Filename ="Ener"
            case (5)
               N = Latt%N*Latt_unit%N_coord;   Filename = "Kin_b" 
           case default
              Write(6,*) ' Error in Alloc_obs '
           end select
           Call Obser_Vec_make(Obs_scal(I),N,Filename)
        enddo
        
        ! Equal time correlators
        Allocate ( Obs_eq(9) )
        Do I = 1,Size(Obs_eq,1)
           select case (I)
           case (1)
              Filename = "Green"
           case (2)
              Filename = "SpinZ"
           case (3)
              Filename = "SpinT"
           case (4)
              Filename = "Den"
           case (5)
              Filename = "QSH"
           case (6)
              Filename = "SC"
           case (7)
              Filename = "SCSUN"
           case (8)
              Filename = "Kin"
           case (9)
              Filename = "SU2col"
           case default
              Write(6,*) ' Error in Alloc_obs '
           end select
           Nt = 1
           Channel = '--'
           if ( I == 5 ) then
              Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_Hex, Channel, dtau)
           elseif ( I == 8) then
              Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_Unit_nnBond, Channel, dtau)
           else
              Call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_Obs, Channel, dtau)
           endif
        enddo


        If (Ltau == 1) then
           ! Time-displaced correlators
           Allocate ( Obs_tau(9) )
           Do I = 1,Size(Obs_tau,1)
              select case (I)
              case (1)
                 Channel = 'P' ; Filename = "Green"
              case (2)
                 Channel = 'PH'; Filename = "SpinZ"
              case (3)
                 Channel = 'PH'; Filename = "SpinT"
              case (4)
                 Channel = 'PH'; Filename = "Den"
              case (5)
                 Channel = 'PH'; Filename = "QSH"
              case (6)
                 Channel = 'PH'; Filename = "SC"
              case (7)
                 Channel = 'PH'; Filename = "SCSUN"
              case (8)
                 Channel = 'PH'; Filename = "Kin"
              case (9)
                 Channel = 'PH'; Filename = "SU2col"
              case default
                 Write(6,*) ' Error in Alloc_obs '
              end select
              Nt = Ltrot+1-2*Thtrot
              If(Projector) Channel = 'T0'
              if (I == 5 ) then
                 Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_Hex, Channel, dtau)
              elseif ( I == 8) then
                 Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_Unit_nnBond, Channel, dtau)
              else
                 Call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_Obs, Channel, dtau)
              endif
           enddo
        endif

      End Subroutine Alloc_obs

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes equal time observables
!> @details
!> @param [IN] Gr   Complex(:,:,:)
!> \verbatim
!>  Green function: Gr(I,J,nf) = <c_{I,nf } c^{dagger}_{J,nf } > on time slice ntau
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!> @param [IN] Ntau Integer
!> \verbatim
!>  Time slice
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

        Use Predefined_Obs
 
        Implicit none
 
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
        Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
        Integer, INTENT(IN)          :: Ntau
        Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
 
        !Local
        Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZK
        Complex (Kind=Kind(0.d0)) :: Z, Zrho, Zkin, ZPot, ZSC, Zn, ZP,ZS, ZZ, ZX, ZY, ZDen, Zcol
        Integer :: I,J, imj, nf,  I1, I2, J1, J2, ns, no, no1, no2, nc, I_up, I_dn, J_up, J_dn, no_i, no_j, norb_I, norb_J, I1_u, I1_d, J1_u, J1_d, I2_u,I2_d,J2_u,J2_d,nb_I,nb_J
        Real    (Kind=Kind(0.d0)) :: X
        Real(8)  ::   X_I, X_J
 
        ZP = PHASE/Real(Phase, kind(0.D0))
        ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
 
        ZS = ZS*Mc_step_weight
                   
         
 
        Do nf = 1,N_FL
           Do I = 1,Ndim
              Do J = 1,Ndim
                 GRC(I, J, nf) = -GR(J, I, nf)
              Enddo
              GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
           Enddo
        Enddo
        ! GRC(i,j,nf) = < c^{dagger}_{i,nf } c_{j,nf } >
 
        ! Compute scalar observables.
        Do I = 1,Size(Obs_scal,1)
           Obs_scal(I)%N         =  Obs_scal(I)%N + 1
           Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
        Enddo
 
 
        Zkin = cmplx(0.d0, 0.d0, kind(0.D0))
        Call Predefined_Hoppings_Compute_Kin(Hopping_Matrix,List,Invlist, Latt, Latt_unit, GRC, ZKin)
        Zkin = Zkin* dble(N_SUN)
        Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS
 
 
        Zn=cmplx(dble(N_sun),0.d0,kind(0.d0))
        ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
        Do nf = 1,N_FL
           Do I = 1,Latt%N
              Do no1 = 1, 6  ! Sum over bonds in hexagon 
                 Do no2 = 1, 6
                    ! \nu_{ij}\nu_{kl}=1 for all bonds in Invlist_Bonds. Negative bonds covered by h.c.
                    Zpot   = Zpot + Predefined_Obs_SO3_eq(I, I, no1, no2, GR, GRC, Invlist_Bonds, N_SUN, N_FL)
                 enddo
              enddo
           Enddo
        Enddo
        Zpot=-Zpot*cmplx(Ham_lambda,0.d0,Kind(0.d0)) ! *Zn (due to color symmetry) /Zn (in Hamiltonian) -> cancels
        Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


        Zrho = cmplx(0.d0,0.d0, kind(0.D0))
        Do nf = 1,N_FL
           Do I = 1,Ndim
              Zrho = Zrho + GRC(I,I,nf)
           enddo
        enddo
        Zrho = Zrho* dble(N_SUN)
        Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
        Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS

        Zn = cmplx(dble(N_SUN), 0.d0, kind(0.d0))
          nc = 0
          Do i = 1, Latt%N
            do no = 1, 3
               nc = nc + 1
               do ns = 1,2 ! Spin
                
                  I1 = Invlist_nnBonds(I,no,ns)
                  J1 = Invlist_nnBonds(I,no,ns+2)
                
                  do nf = 1, N_FL
                     Obs_scal(5)%Obs_vec(nc)  =   Obs_scal(5)%Obs_vec(nc) + Zn*(GRC(I1,J1,nf) + GRC(J1,I1,nf))*ZP*ZS
                  enddo
               enddo
            enddo
          enddo

        
        ! Standard two-point correlations
        ! Call Predefined_Obs_eq_Green_measure  ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(1) )
        ! Call Predefined_Obs_eq_SpinSUN_measure( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(2) )
        ! Call Predefined_Obs_eq_Den_measure    ( Latt, Latt_unit, List,  GR, GRC, N_SUN, ZS, ZP, Obs_eq(3) )


        Do I = 1,Size(Obs_eq,1)
           Obs_eq(I)%N         =  Obs_eq(I)%N + 1
           Obs_eq(I)%Ave_sign  =  Obs_eq(I)%Ave_sign + Real(ZS,kind(0.d0))
        Enddo

        ! QSH correlations
        Do I = 1,Latt%N
           Do no1  = 1, 6
              Do J = 1,Latt%N
                 Imj = latt%imj(I,J)
                 Do no2  = 1, 6
                    Z   = Predefined_Obs_SO3_eq( I, J, no1, no2, GR, GRC, Invlist_Bonds, N_SUN, N_FL)
                    Obs_eq(5)%Obs_Latt(imj,1,no1,no2) =  Obs_eq(5)%Obs_Latt(imj,1,no1,no2) + Z*ZP*ZS
                 enddo
              enddo
           enddo
        enddo

        ! Kin correlations
        Do I = 1,Latt%N 
           Do no1 = 1,3
              Do J = 1,Latt%N 
                 Imj = latt%imj(I,J)
                 Do no2 = 1,3
                    Z = Predefined_Obs_Kin_eq( I, J, no1, no2, GR, GRC, Invlist_nnBonds, N_SUN, N_FL)
                    Obs_eq(8)%Obs_Latt(imj,1,no1,no2) = Obs_eq(8)%Obs_Latt(imj,1,no1,no2) + Z*ZP*ZS
                 enddo
              enddo
              Z = Predefined_Obs_Kin_eq0( I, no1, GR, GRC, Invlist_nnBonds, N_SUN, N_FL)
              Obs_eq(8)%Obs_Latt0(no1) = Obs_eq(8)%Obs_Latt0(no1) + Z*ZP*ZS
           enddo
        enddo

        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        Do I = 1,Latt%N 
           Do no1 = 1,2
              I1 = Invlist(I,no1)     ! Orbital no1 of unit cell I spin up
              I2 = Invlist(I,no1+2)   ! Orbital no1 of unit cell I spin down
              Do J = 1,Latt%N 
                 imj = latt%imj(I,J)
                 Do no2 = 1,2
                    J1 = Invlist(J,no2)
                    J2 = Invlist(J,no2+2)

                    ZX = (GRC(I1,I2,1) * GRC(J1,J2,1) + GRC(I1,I2,1) * GRC(J2,J1,1) &
                       +  GRC(I2,I1,1) * GRC(J1,J2,1) + GRC(I2,I1,1) * GRC(J2,J1,1))*ZN*ZN &
                       + (GRC(I1,J2,1) * GR(I2,J1,1)  + GRC(I1,J1,1) * GR(I2,J2,1) &
                       +  GRC(I2,J2,1) * GR(I1,J1,1)  + GRC(I2,J1,1) * GR(I1,J2,1))*ZN
                    ZY = (-GRC(I1,I2,1) * GRC(J1,J2,1) + GRC(I1,I2,1) * GRC(J2,J1,1) &
                       +   GRC(I2,I1,1) * GRC(J1,J2,1) - GRC(I2,I1,1) * GRC(J2,J1,1))*ZN*ZN &
                       + (-GRC(I1,J2,1) * GR(I2,J1,1)  + GRC(I1,J1,1) * GR(I2,J2,1) &
                       +   GRC(I2,J2,1) * GR(I1,J1,1)  - GRC(I2,J1,1) * GR(I1,J2,1))*ZN
                    ZZ = (GRC(I1,J1,1) * GR(I1,J1,1) + GRC(I2,J2,1) * GR(I2,J2,1) &
                       -  GRC(I1,J2,1) * GR(I1,J2,1) - GRC(I2,J1,1) * GR(I2,J1,1))*ZN &
                       + (GRC(I2,I2,1) - GRC(I1,I1,1))*(GRC(J2,J2,1) - GRC(J1,J1,1))*ZN*ZN

                    ZDen = (GRC(I1,I1,1) + GRC(I2,I2,1)) * (GRC(J1,J1,1) + GRC(J2,J2,1))*ZN*ZN &
                       + (GRC(I1,J1,1) * GR(I1,J1,1)   +  GRC(I2,J2,1) * GR(I2,J2,1) &
                       + GRC(I1,J2,1) * GR(I1,J2,1) + GRC(I2,J1,1) * GR(I2,J1,1))*ZN

                    ! SC pair-pair correlations

                    ZSC = GRC(I1,J1,1)*GRC(I2,J2,1) + GRC(J1,I1,1)*GRC(J2,I2,1) &
                       - GRC(I1,J2,1)*GRC(I2,J1,1) - GRC(J1,I2,1)*GRC(J2,I1,1)
                    ZSC = ZSC * cmplx(dble(N_SUN)*0.25d0, 0.d0, kind(0.D0))

                    Z = GRC(I1,J1,1)*GRC(I2,J2,1) - GRC(I1,J2,1)*GRC(I2,J1,1) ! SC pair-pair without h.c.

                    Zcol = 6.d0 * ( GRC(I1,J1,1)*GR(I1,J1,1) + GRC(I1,J2,1)*GR(I1,J2,1) &
                                  + GRC(I2,J1,1)*GR(I2,J1,1) + GRC(I2,J2,1)*GR(I2,J2,1) )

                    Obs_eq(1)%Obs_Latt(imj,1,no1,no2) = Obs_eq(1)%Obs_Latt(imj,1,no1,no2) + (GRC(I1,J1,1) + GRC(I2,J2,1))*ZN*ZP*ZS
                    Obs_eq(2)%Obs_Latt(imj,1,no1,no2) = Obs_eq(2)%Obs_Latt(imj,1,no1,no2) +  ZZ *ZP*ZS
                    Obs_eq(3)%Obs_Latt(imj,1,no1,no2) = Obs_eq(3)%Obs_Latt(imj,1,no1,no2) + (ZX + ZY + ZZ)*ZP*ZS/3.d0
                    Obs_eq(4)%Obs_Latt(imj,1,no1,no2) = Obs_eq(4)%Obs_Latt(imj,1,no1,no2) +  ZDen * ZP * ZS
                    Obs_eq(6)%Obs_Latt(imj,1,no1,no2) = Obs_eq(6)%Obs_Latt(imj,1,no1,no2) + ZSC*ZP*ZS
                    Obs_eq(7)%Obs_Latt(imj,1,no1,no2) = Obs_eq(7)%Obs_Latt(imj,1,no1,no2) + (Z**ZN)*ZP*ZS
                    Obs_eq(9)%Obs_Latt(imj,1,no1,no2) = Obs_eq(9)%Obs_Latt(imj,1,no1,no2) + (Zcol)*ZP*ZS
                 enddo
              enddo
              Obs_eq(4)%Obs_Latt0(no1) = Obs_eq(4)%Obs_Latt0(no1) + (GRC(I1,I1,1) + GRC(I2,I2,1))*ZN *  ZP*ZS
           enddo
        enddo

      end Subroutine Obser
!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Computes time displaced  observables
!> @details
!> @param [IN] NT, Integer
!> \verbatim
!>  Imaginary time
!> \endverbatim
!> @param [IN] GT0, GTT, G00, GTT,  Complex(:,:,:)
!> \verbatim
!>  Green functions:
!>  GT0(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(0  )>
!>  G0T(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(tau)>
!>  G00(I,J,nf) = <T c_{I,nf }(0  ) c^{dagger}_{J,nf }(0  )>
!>  GTT(I,J,nf) = <T c_{I,nf }(tau) c^{dagger}_{J,nf }(tau)>
!> \endverbatim
!> @param [IN] Phase   Complex
!> \verbatim
!>  Phase
!> \endverbatim
!-------------------------------------------------------------------
      Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE, Mc_step_weight )

        Use Predefined_Obs

        Implicit none

        Integer         , INTENT(IN) :: NT
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
        Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
        Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

        !Locals
        Complex (Kind=Kind(0.d0)) :: Z, Zn, ZP, ZS, ZSC, ZZ, ZX, ZY, ZDen, Z1,Z2,Z3,Z4,Zcol
        Real    (Kind=Kind(0.d0)) :: X
        Real(8)  ::   X_I, X_J
        Integer :: IMJ, I, J, I1, J1, I2, J2, no1, no2, I_up, I_dn, J_up, J_dn, no_I, no_j, norb_I, norb_J, Sz, nb_I,nb_J,I1_u,J1_u,I1_d,J1_d,I2_u,J2_u,I2_d,J2_d

        ZP = PHASE/Real(Phase, kind(0.D0))
        ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

        ZS = ZS*Mc_step_weight

         
        ! Standard two-point correlations

        ! Call Predefined_Obs_tau_Green_measure  ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(1) )
        ! Call Predefined_Obs_tau_SpinSUN_measure( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(2) )
        ! Call Predefined_Obs_tau_Den_measure    ( Latt, Latt_unit, List, NT, GT0,G0T,G00,GTT,  N_SUN, ZS, ZP, Obs_tau(3) )

        
        If (NT == 0 ) then
           Do I = 1,Size(Obs_tau,1)
              Obs_tau(I)%N         =  Obs_tau(I)%N + 1
              Obs_tau(I)%Ave_sign  =  Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
           Enddo
        Endif
        
        ! QSH correlations
        Do I = 1,Latt%N
           Do no1  = 1, 6
              Do J = 1,Latt%N
                 Imj = latt%imj(I,J)
                 Do no2  = 1, 6
                    Z   = Predefined_Obs_SO3_tau( I, J, no1, no2, GT0,G0T,G00,GTT, Invlist_Bonds, N_SUN, N_FL)
                    Obs_tau(5)%Obs_Latt(imj,NT+1,no1,no2) =  Obs_tau(5)%Obs_Latt(imj,NT+1,no1,no2) + Z*ZP*ZS
                 enddo
              enddo
           enddo
        enddo

        ! Kin correlations
        Do I = 1,Latt%N 
           Do no1 = 1,3
              Do J = 1,Latt%N 
                 Imj = latt%imj(I,J)
                 Do no2 = 1,3
                    Z = Predefined_Obs_Kin_tau( I, J, no1, no2, GT0, G0T, G00, GTT, Invlist_nnBonds, N_SUN, N_FL)
                    Obs_tau(8)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(8)%Obs_Latt(imj,NT+1,no1,no2) + Z*ZP*ZS
                 enddo
              enddo
              Z = Predefined_Obs_Kin_tau0( I, no1, GT0, G0T, G00, GTT, Invlist_nnBonds, N_SUN, N_FL)
              Obs_tau(8)%Obs_Latt0(no1) = Obs_tau(8)%Obs_Latt0(no1) + Z*ZP*ZS
           enddo
        enddo

        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        Do I = 1,Latt%N 
           Do no1 = 1,2
              I1 = Invlist(I,no1)     ! Orbital no1 of unit cell I spin up
              I2 = Invlist(I,no1+2)   ! Orbital no1 of unit cell I spin down
              Do J = 1,Latt%N 
                 imj = latt%imj(I,J)
                 Do no2 = 1,2
                    J1 = Invlist(J,no2)
                    J2 = Invlist(J,no2+2)
                    
                    ZZ = (GTT(I1,I1,1) -  GTT(I2,I2,1) ) * ( G00(J1,J1,1)  -  G00(J2,J2,1) ) *ZN*ZN   &
                       +(-G0T(J1,I1,1) * GT0(I1,J1,1)  -  G0T(J2,I2,1) * GT0(I2,J2,1) &
                       +  G0T(J2,I1,1) * GT0(I1,J2,1)  +  G0T(J1,I2,1) * GT0(I2,J1,1))*ZN
                    ZX = (GTT(I2,I1,1) * G00(J2,J1,1) + GTT(I2,I1,1) * G00(J1,J2,1) &
                       +  GTT(I1,I2,1) * G00(J2,J1,1) + GTT(I1,I2,1) * G00(J1,J2,1))*ZN*ZN &
                       +(-G0T(J2,I1,1) * GT0(I2,J1,1) - G0T(J1,I1,1) * GT0(I2,J2,1) &
                       -  G0T(J2,I2,1) * GT0(I1,J1,1) - G0T(J1,I2,1) * GT0(I1,J2,1))*ZN
                    ZY =(-GTT(I2,I1,1) * G00(J2,J1,1) + GTT(I2,I1,1) * G00(J1,J2,1) &
                       +  GTT(I1,I2,1) * G00(J2,J1,1) - GTT(I1,I2,1) * G00(J1,J2,1))*ZN*ZN &
                       + (G0T(J2,I1,1) * GT0(I2,J1,1) - G0T(J1,I1,1) * GT0(I2,J2,1) &
                       -  G0T(J2,I2,1) * GT0(I1,J1,1) + G0T(J1,I2,1) * GT0(I1,J2,1))*ZN

                    ZDen = (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I1,I1,1) - GTT(I2,I2,1) ) * &
                       &   (cmplx(2.d0,0.d0,kind(0.d0)) -  G00(J1,J1,1) - G00(J2,J2,1) )*ZN*ZN   &
                       & + (-G0T(J1,I1,1) * GT0(I1,J1,1) - G0T(J2,I2,1) * GT0(I2,J2,1) &
                       &    -G0T(J2,I1,1) * GT0(I1,J2,1) - G0T(J1,I2,1) * GT0(I2,J1,1))*ZN

                    ! SC pair-pair correlations
                    ZSC = G0T(J1,I1,1)*G0T(J2,I2,1) - G0T(J1,I2,1)*G0T(J2,I1,1) &
                       & + GT0(I2,J2,1)*GT0(I1,J1,1) - GT0(I2,J1,1)*GT0(I1,J2,1)
                    ZSC = ZSC * cmplx(dble(N_SUN)*0.25d0, 0.d0, kind(0.D0))

                    Z = G0T(J1,I1,1)*G0T(J2,I2,1) - G0T(J1,I2,1)*G0T(J2,I1,1)

                    Zcol = - 6.d0 * ( G0T(J1,I1,1)*GT0(I1,J1,1) + G0T(J2,I1,1)*GT0(I1,J2,1) &
                                  + G0T(J1,I2,1)*GT0(I2,J1,1) + G0T(J2,I1,1)*GT0(I2,J2,1) )
                    
                    Obs_tau(1)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(1)%Obs_Latt(imj,NT+1,no1,no2) + (GT0(I1,J1,1) + GT0(I2,J2,1))*ZN*ZP*ZS
                    Obs_tau(2)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(2)%Obs_Latt(imj,NT+1,no1,no2) +  ZZ  *ZP*ZS
                    Obs_tau(3)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(3)%Obs_Latt(imj,NT+1,no1,no2) + (Zx + ZY + ZZ)*ZP*ZS/3.d0
                    Obs_tau(4)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(4)%Obs_Latt(imj,NT+1,no1,no2) +  ZDen * ZP * ZS

                    Obs_tau(6)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(6)%Obs_Latt(imj,NT+1,no1,no2) + ZSC*ZP*ZS
                    Obs_tau(7)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(7)%Obs_Latt(imj,NT+1,no1,no2) + (Z**ZN)*ZP*ZS
                    Obs_tau(9)%Obs_Latt(imj,NT+1,no1,no2) = Obs_tau(9)%Obs_Latt(imj,NT+1,no1,no2) + (Zcol)*ZP*ZS
                 enddo
              enddo
              Obs_tau(4)%Obs_Latt0(no1) = Obs_tau(4)%Obs_Latt0(no1) + &
                    &  (cmplx(2.d0,0.d0,kind(0.d0)) -  GTT(I1,I1,1) - GTT(I2,I2,1))*ZN  *  ZP*ZS

           enddo
        enddo

      end Subroutine OBSERT 


      !-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure SO(3) spin-spin.
!>  Let O_{r,\delta} = ic^{\dagger}_{i}\vec{\sigma}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <O_{r,\delta} O_{r^{\prime},\delta^{\prime}} >
!-----------------------------------------------------------------------
     Complex (Kind=Kind(0.d0)) function Predefined_Obs_SO3_eq(XR, XRP, nb1, nb2, GR, GRC, Invlist_Bonds, N_SUN, N_FL)
  
        Implicit none
        
        Integer, Intent(IN) ::  XR, XRP, nb1, nb2, N_SUN, N_FL
        Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_Bonds
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GR,GRC
     
        ! local
        Integer ::  I0_up,J0_up,K0_up,L0_up, I0_dn,J0_dn,K0_dn,L0_dn, N, Isl
        Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Z2, Z3, Z4, Zsgn
        Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2), GC_kl(2,2), GC_il(2,2), G_jk(2,2)

        Z = cmplx(0.d0,0.d0,Kind(0.d0))
        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        if (N_FL == 1 ) then
           do Isl = 1, 4
              select case (Isl)
              case (1)
                    I0_up = Invlist_Bonds(XR , nb1, 1); I0_dn = Invlist_Bonds(XR , nb1, 2)
                    J0_up = Invlist_Bonds(XR , nb1, 3); J0_dn = Invlist_Bonds(XR , nb1, 4)
                    k0_up = Invlist_Bonds(XRP, nb2, 1); k0_dn = Invlist_Bonds(XRP, nb2, 2)
                    l0_up = Invlist_Bonds(XRP, nb2, 3); l0_dn = Invlist_Bonds(XRP, nb2, 4)
                    Zsgn = cmplx(-1.d0,0.d0,Kind(0.d0))
              case (2)
                    I0_up = Invlist_Bonds(XR , nb1, 1); I0_dn = Invlist_Bonds(XR , nb1, 2)
                    J0_up = Invlist_Bonds(XR , nb1, 3); J0_dn = Invlist_Bonds(XR , nb1, 4)
                    k0_up = Invlist_Bonds(XRP, nb2, 3); k0_dn = Invlist_Bonds(XRP, nb2, 4)
                    l0_up = Invlist_Bonds(XRP, nb2, 1); l0_dn = Invlist_Bonds(XRP, nb2, 2)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (3)
                    I0_up = Invlist_Bonds(XR , nb1, 3); I0_dn = Invlist_Bonds(XR , nb1, 4)
                    J0_up = Invlist_Bonds(XR , nb1, 1); J0_dn = Invlist_Bonds(XR , nb1, 2)
                    k0_up = Invlist_Bonds(XRP, nb2, 1); k0_dn = Invlist_Bonds(XRP, nb2, 2)
                    l0_up = Invlist_Bonds(XRP, nb2, 3); l0_dn = Invlist_Bonds(XRP, nb2, 4)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (4)
                    I0_up = Invlist_Bonds(XR , nb1, 3); I0_dn = Invlist_Bonds(XR , nb1, 4)
                    J0_up = Invlist_Bonds(XR , nb1, 1); J0_dn = Invlist_Bonds(XR , nb1, 2)
                    k0_up = Invlist_Bonds(XRP, nb2, 3); k0_dn = Invlist_Bonds(XRP, nb2, 4)
                    l0_up = Invlist_Bonds(XRP, nb2, 1); l0_dn = Invlist_Bonds(XRP, nb2, 2)
                    Zsgn = cmplx(-1.d0,0.d0,Kind(0.d0))
              case default
                 Write(6,*) ' Error in SO3 obs '
              end select
              GC_ij(1,1) = GRC(I0_up,J0_up,1); GC_ij(1,2) = GRC(I0_up,J0_dn,1)
              GC_ij(2,1) = GRC(I0_dn,J0_up,1); GC_ij(2,2) = GRC(I0_dn,J0_dn,1)

              GC_kl(1,1) = GRC(k0_up,l0_up,1); GC_kl(1,2) = GRC(k0_up,l0_dn,1)
              GC_kl(2,1) = GRC(k0_dn,l0_up,1); GC_kl(2,2) = GRC(k0_dn,l0_dn,1)

              GC_il(1,1) = GRC(I0_up,l0_up,1); GC_il(1,2) = GRC(I0_up,l0_dn,1)
              GC_il(2,1) = GRC(I0_dn,l0_up,1); GC_il(2,2) = GRC(I0_dn,l0_dn,1)

              G_jk(1,1)  = GR(j0_up,k0_up,1) ; G_jk(1,2)  = GR(j0_up,k0_dn,1)
              G_jk(2,1)  = GR(j0_dn,k0_up,1) ; G_jk(2,2)  = GR(j0_dn,k0_dn,1)

              Z1 = (2.d0*GC_ij(2,1)*GC_kl(1,2) + 2.d0*GC_ij(1,2)*GC_kl(2,1) + &
                    (GC_ij(1,1) - GC_ij(2,2))*(GC_kl(1,1) - GC_kl(2,2)))*ZN*ZN + &
                    ((GC_il(1,1) + 2.d0*GC_il(2,2))*G_jk(1,1)- GC_il(1,2)*G_jk(1,2) - &
                    GC_il(2,1)*G_jk(2,1) + (2.d0*GC_il(1,1) + GC_il(2,2))*G_jk(2,2))*ZN

              Z  = Z + Z1*Zsgn
           enddo
        endif

        Predefined_Obs_SO3_eq = Z
        
     end function Predefined_Obs_SO3_eq

     !-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure SO(3) spin-spin.
!>  Let O_{r,\delta} = ic^{\dagger}_{i}\vec{\sigma}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <O_{r,\delta} O_{r^{\prime},\delta^{\prime}} >
!-----------------------------------------------------------------------
     Complex (Kind=Kind(0.d0)) function Predefined_Obs_SO3_tau(XR, XRP, nb1, nb2,\
         GT0,G0T,G00,GTT, Invlist_Bonds, N_SUN, N_FL)
 
       Implicit none
       
       Integer, Intent(IN) ::  XR, XRP, nb1, nb2, N_SUN, N_FL
       Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_Bonds
       Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GT0,G0T,G00,GTT
      
       ! local
       Integer ::  I0_up,J0_up,K0_up,L0_up, I0_dn,J0_dn,K0_dn,L0_dn, N, Isl
       Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Z2, Z3, Z4, Zsgn, Delta_1, Delta_2
       Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2), GC_kl(2,2), GC_il(2,2), G_jk(2,2)

       Z = cmplx(0.d0,0.d0,Kind(0.d0))
       ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
       if (N_FL == 1 ) then
           do Isl = 1, 4
               select case (Isl)
               case (1)
                   I0_up = Invlist_Bonds(XR , nb1, 1); I0_dn = Invlist_Bonds(XR , nb1, 2)
                   J0_up = Invlist_Bonds(XR , nb1, 3); J0_dn = Invlist_Bonds(XR , nb1, 4)
                   k0_up = Invlist_Bonds(XRP, nb2, 1); k0_dn = Invlist_Bonds(XRP, nb2, 2)
                   l0_up = Invlist_Bonds(XRP, nb2, 3); l0_dn = Invlist_Bonds(XRP, nb2, 4)
                   Zsgn = cmplx(-1.d0,0.d0,Kind(0.d0))
               case (2)
                   I0_up = Invlist_Bonds(XR , nb1, 1); I0_dn = Invlist_Bonds(XR , nb1, 2)
                   J0_up = Invlist_Bonds(XR , nb1, 3); J0_dn = Invlist_Bonds(XR , nb1, 4)
                   k0_up = Invlist_Bonds(XRP, nb2, 3); k0_dn = Invlist_Bonds(XRP, nb2, 4)
                   l0_up = Invlist_Bonds(XRP, nb2, 1); l0_dn = Invlist_Bonds(XRP, nb2, 2)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case (3)
                   I0_up = Invlist_Bonds(XR , nb1, 3); I0_dn = Invlist_Bonds(XR , nb1, 4)
                   J0_up = Invlist_Bonds(XR , nb1, 1); J0_dn = Invlist_Bonds(XR , nb1, 2)
                   k0_up = Invlist_Bonds(XRP, nb2, 1); k0_dn = Invlist_Bonds(XRP, nb2, 2)
                   l0_up = Invlist_Bonds(XRP, nb2, 3); l0_dn = Invlist_Bonds(XRP, nb2, 4)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case (4)
                   I0_up = Invlist_Bonds(XR , nb1, 3); I0_dn = Invlist_Bonds(XR , nb1, 4)
                   J0_up = Invlist_Bonds(XR , nb1, 1); J0_dn = Invlist_Bonds(XR , nb1, 2)
                   k0_up = Invlist_Bonds(XRP, nb2, 3); k0_dn = Invlist_Bonds(XRP, nb2, 4)
                   l0_up = Invlist_Bonds(XRP, nb2, 1); l0_dn = Invlist_Bonds(XRP, nb2, 2)
                   Zsgn = cmplx(-1.d0,0.d0,Kind(0.d0))
               case default
                  Write(6,*) ' Error in SO3 obs '
               end select

               Delta_1 = cmplx(0.d0,0.d0,Kind(0.d0))
               Delta_2 = cmplx(0.d0,0.d0,Kind(0.d0))
               if (I0_up .eq. J0_up) Delta_1 = cmplx(1.d0,0.d0,Kind(0.d0))
               if (k0_up .eq. l0_up) Delta_2 = cmplx(1.d0,0.d0,Kind(0.d0))

               GC_ij(1,1) = Delta_1-GTT(J0_up,I0_up,1); GC_ij(1,2) = -GTT(J0_dn,I0_up,1)
               GC_ij(2,1) = -GTT(J0_up,I0_dn,1); GC_ij(2,2) = Delta_1-GTT(J0_dn,I0_dn,1)

               GC_kl(1,1) = Delta_2-G00(l0_up,k0_up,1); GC_kl(1,2) = -G00(l0_dn,k0_up,1)
               GC_kl(2,1) = -G00(l0_up,k0_dn,1); GC_kl(2,2) = Delta_2-G00(l0_dn,k0_dn,1)

               GC_il(1,1) = -G0T(l0_up,i0_up,1); GC_il(1,2) = -G0T(l0_dn,i0_up,1)
               GC_il(2,1) = -G0T(l0_up,i0_dn,1); GC_il(2,2) = -G0T(l0_dn,i0_dn,1)

               G_jk(1,1)  = GT0(j0_up,k0_up,1) ; G_jk(1,2)  = GT0(j0_up,k0_dn,1)
               G_jk(2,1)  = GT0(j0_dn,k0_up,1) ; G_jk(2,2)  = GT0(j0_dn,k0_dn,1)

               Z1 = (2.d0*GC_ij(2,1)*GC_kl(1,2) + 2.d0*GC_ij(1,2)*GC_kl(2,1) + &
                   (GC_ij(1,1) - GC_ij(2,2))*(GC_kl(1,1) - GC_kl(2,2)))*ZN*ZN + &
                   ((GC_il(1,1) + 2.d0*GC_il(2,2))*G_jk(1,1)- GC_il(1,2)*G_jk(1,2) - &
                   GC_il(2,1)*G_jk(2,1) + (2.d0*GC_il(1,1) + GC_il(2,2))*G_jk(2,2))*ZN

               Z  = Z + Z1*Zsgn
           enddo
       endif

       Predefined_Obs_SO3_tau = Z
       
     end function Predefined_Obs_SO3_tau

!-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure Kin-Kin.
!>  Let K_{r,\delta} = ic^{\dagger}_{i}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <K_{r,\delta} K_{r^{\prime},\delta^{\prime}} >
!-----------------------------------------------------------------------
     Complex (Kind=Kind(0.d0)) function Predefined_Obs_Kin_eq(XR, XRP, nb1, nb2, GR, GRC, Invlist_nnBonds, N_SUN, N_FL)
  
        Implicit none
        
        Integer, Intent(IN) ::  XR, XRP, nb1, nb2, N_SUN, N_FL
        Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_nnBonds
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GR,GRC
     
        ! local
        Integer ::  I0_up,J0_up,K0_up,L0_up, I0_dn,J0_dn,K0_dn,L0_dn, N, Isl, ns1, ns2
        Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Zsgn
        Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2), GC_kl(2,2), GC_il(2,2), G_jk(2,2)

        Z = cmplx(0.d0,0.d0,Kind(0.d0))
        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        if (N_FL == 1 ) then
           do Isl = 1, 4
              select case (Isl)
              case (1)
                    I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                    J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                    k0_up = Invlist_nnBonds(XRP, nb2, 1); k0_dn = Invlist_nnBonds(XRP, nb2, 2)
                    l0_up = Invlist_nnBonds(XRP, nb2, 3); l0_dn = Invlist_nnBonds(XRP, nb2, 4)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (2)
                    I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                    J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                    k0_up = Invlist_nnBonds(XRP, nb2, 3); k0_dn = Invlist_nnBonds(XRP, nb2, 4)
                    l0_up = Invlist_nnBonds(XRP, nb2, 1); l0_dn = Invlist_nnBonds(XRP, nb2, 2)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (3)
                    I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                    J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                    k0_up = Invlist_nnBonds(XRP, nb2, 1); k0_dn = Invlist_nnBonds(XRP, nb2, 2)
                    l0_up = Invlist_nnBonds(XRP, nb2, 3); l0_dn = Invlist_nnBonds(XRP, nb2, 4)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (4)
                    I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                    J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                    k0_up = Invlist_nnBonds(XRP, nb2, 3); k0_dn = Invlist_nnBonds(XRP, nb2, 4)
                    l0_up = Invlist_nnBonds(XRP, nb2, 1); l0_dn = Invlist_nnBonds(XRP, nb2, 2)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case default
                 Write(6,*) ' Error in Kin obs '
              end select
              GC_ij(1,1) = GRC(I0_up,J0_up,1); GC_ij(1,2) = GRC(I0_up,J0_dn,1)
              GC_ij(2,1) = GRC(I0_dn,J0_up,1); GC_ij(2,2) = GRC(I0_dn,J0_dn,1)

              GC_kl(1,1) = GRC(k0_up,l0_up,1); GC_kl(1,2) = GRC(k0_up,l0_dn,1)
              GC_kl(2,1) = GRC(k0_dn,l0_up,1); GC_kl(2,2) = GRC(k0_dn,l0_dn,1)

              GC_il(1,1) = GRC(I0_up,l0_up,1); GC_il(1,2) = GRC(I0_up,l0_dn,1)
              GC_il(2,1) = GRC(I0_dn,l0_up,1); GC_il(2,2) = GRC(I0_dn,l0_dn,1)

              G_jk(1,1)  = GR(j0_up,k0_up,1) ; G_jk(1,2)  = GR(j0_up,k0_dn,1)
              G_jk(2,1)  = GR(j0_dn,k0_up,1) ; G_jk(2,2)  = GR(j0_dn,k0_dn,1)
              do ns1 = 1,2
                 do ns2 = 1,2
                    Z1 = (GC_ij(ns1,ns1)*GC_kl(ns2,ns2))*ZN*ZN + (GC_il(ns1,ns2)*G_jk(ns1,ns2))*ZN

                    Z  = Z + Z1*Zsgn
                 enddo
              enddo
           enddo
        endif

        Predefined_Obs_Kin_eq = Z
        
     end function Predefined_Obs_Kin_eq

!-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure Kin-background.
!>  Let K_{r,\delta} = ic^{\dagger}_{i}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <K_{r,\delta}>
!-----------------------------------------------------------------------

     Complex (Kind=Kind(0.d0)) function Predefined_Obs_Kin_eq0(XR, nb1, GR, GRC, Invlist_nnBonds, N_SUN, N_FL)
  
        Implicit none
        
        Integer, Intent(IN) ::  XR, nb1, N_SUN, N_FL
        Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_nnBonds
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GR,GRC
     
        ! local
        Integer ::  I0_up,J0_up, I0_dn,J0_dn, N, Isl, ns1
        Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Zsgn
        Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2)

        Z = cmplx(0.d0,0.d0,Kind(0.d0))
        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        if (N_FL == 1 ) then
           do Isl = 1, 2
              select case (Isl)
              case (1)
                    I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                    J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (2)
                    I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                    J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case default
                 Write(6,*) ' Error in Kin obs '
              end select
              GC_ij(1,1) = GRC(I0_up,J0_up,1); GC_ij(1,2) = GRC(I0_up,J0_dn,1)
              GC_ij(2,1) = GRC(I0_dn,J0_up,1); GC_ij(2,2) = GRC(I0_dn,J0_dn,1)

              do ns1 = 1,2
                 Z1 = GC_ij(ns1,ns1)*ZN

                 Z  = Z + Z1*Zsgn
              enddo
           enddo
        endif

        Predefined_Obs_Kin_eq0 = Z
        
     end function Predefined_Obs_Kin_eq0

!-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure Kin-Kin.
!>  Let K_{r,\delta} = ic^{\dagger}_{i}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <K_{r,\delta} K_{r^{\prime},\delta^{\prime}} >
!-----------------------------------------------------------------------
     Complex (Kind=Kind(0.d0)) function Predefined_Obs_Kin_tau(XR, XRP, nb1, nb2,\
         GT0,G0T,G00,GTT, Invlist_nnBonds, N_SUN, N_FL)
 
       Implicit none
       
       Integer, Intent(IN) ::  XR, XRP, nb1, nb2, N_SUN, N_FL
       Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_nnBonds
       Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GT0,G0T,G00,GTT
      
       ! local
       Integer ::  I0_up,J0_up,K0_up,L0_up, I0_dn,J0_dn,K0_dn,L0_dn, N, Isl, ns1, ns2
       Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Zsgn, Delta_1, Delta_2
       Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2), GC_kl(2,2), GC_il(2,2), G_jk(2,2)

       Z = cmplx(0.d0,0.d0,Kind(0.d0))
       ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
       if (N_FL == 1 ) then
           do Isl = 1, 4
               select case (Isl)
               case (1)
                   I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                   J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                   k0_up = Invlist_nnBonds(XRP, nb2, 1); k0_dn = Invlist_nnBonds(XRP, nb2, 2)
                   l0_up = Invlist_nnBonds(XRP, nb2, 3); l0_dn = Invlist_nnBonds(XRP, nb2, 4)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case (2)
                   I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                   J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                   k0_up = Invlist_nnBonds(XRP, nb2, 3); k0_dn = Invlist_nnBonds(XRP, nb2, 4)
                   l0_up = Invlist_nnBonds(XRP, nb2, 1); l0_dn = Invlist_nnBonds(XRP, nb2, 2)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case (3)
                   I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                   J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                   k0_up = Invlist_nnBonds(XRP, nb2, 1); k0_dn = Invlist_nnBonds(XRP, nb2, 2)
                   l0_up = Invlist_nnBonds(XRP, nb2, 3); l0_dn = Invlist_nnBonds(XRP, nb2, 4)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case (4)
                   I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                   J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                   k0_up = Invlist_nnBonds(XRP, nb2, 3); k0_dn = Invlist_nnBonds(XRP, nb2, 4)
                   l0_up = Invlist_nnBonds(XRP, nb2, 1); l0_dn = Invlist_nnBonds(XRP, nb2, 2)
                   Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
               case default
                  Write(6,*) ' Error in Kin obs '
               end select

               Delta_1 = cmplx(0.d0,0.d0,Kind(0.d0))
               Delta_2 = cmplx(0.d0,0.d0,Kind(0.d0))
               if (I0_up .eq. J0_up) Delta_1 = cmplx(1.d0,0.d0,Kind(0.d0))
               if (k0_up .eq. l0_up) Delta_2 = cmplx(1.d0,0.d0,Kind(0.d0))

               GC_ij(1,1) = Delta_1-GTT(J0_up,I0_up,1); GC_ij(1,2) = -GTT(J0_dn,I0_up,1)
               GC_ij(2,1) = -GTT(J0_up,I0_dn,1); GC_ij(2,2) = Delta_1-GTT(J0_dn,I0_dn,1)

               GC_kl(1,1) = Delta_2-G00(l0_up,k0_up,1); GC_kl(1,2) = -G00(l0_dn,k0_up,1)
               GC_kl(2,1) = -G00(l0_up,k0_dn,1); GC_kl(2,2) = Delta_2-G00(l0_dn,k0_dn,1)

               GC_il(1,1) = -G0T(l0_up,i0_up,1); GC_il(1,2) = -G0T(l0_dn,i0_up,1)
               GC_il(2,1) = -G0T(l0_up,i0_dn,1); GC_il(2,2) = -G0T(l0_dn,i0_dn,1)

               G_jk(1,1)  = GT0(j0_up,k0_up,1) ; G_jk(1,2)  = GT0(j0_up,k0_dn,1)
               G_jk(2,1)  = GT0(j0_dn,k0_up,1) ; G_jk(2,2)  = GT0(j0_dn,k0_dn,1)

              do ns1 = 1,2
                 do ns2 = 1,2
                    Z1 = (GC_ij(ns1,ns1)*GC_kl(ns2,ns2))*ZN*ZN + (GC_il(ns1,ns2)*G_jk(ns1,ns2))*ZN

                    Z  = Z + Z1*Zsgn
                 enddo
              enddo
           enddo
       endif

       Predefined_Obs_Kin_tau = Z
       
     end function Predefined_Obs_Kin_tau

!-----------------------------------------------------------------------
!> @author
!> ALF-project
!
!>  @brief
!>  Measure Kin-background.
!>  Let K_{r,\delta} = ic^{\dagger}_{i}c_{i+\delta} + h.c.
!>  Routine returns:
!>  <K_{r,\delta}>
!-----------------------------------------------------------------------

     Complex (Kind=Kind(0.d0)) function Predefined_Obs_Kin_tau0(XR, nb1, GT0,G0T,G00,GTT, Invlist_nnBonds, N_SUN, N_FL)
  
        Implicit none
        
        Integer, Intent(IN) ::  XR, nb1, N_SUN, N_FL
        Integer, Dimension(:,:,:), INTENT(IN) :: Invlist_nnBonds
        Complex (Kind=Kind(0.d0)), Dimension(:,:,:), INTENT(IN) :: GT0,G0T,G00,GTT
     
        ! local
        Integer ::  I0_up,J0_up, I0_dn,J0_dn, N, Isl, ns1
        Complex (Kind=Kind(0.d0))  ::  Z, ZK, ZN, Z1, Zsgn, Delta_1
        Complex (Kind=Kind(0.d0))  ::  GC_ij(2,2)

        Z = cmplx(0.d0,0.d0,Kind(0.d0))
        ZN = cmplx(Real(N_SUN,Kind(0.d0)),0.d0,Kind(0.d0))
        if (N_FL == 1 ) then
           do Isl = 1, 2
              select case (Isl)
              case (1)
                    I0_up = Invlist_nnBonds(XR , nb1, 1); I0_dn = Invlist_nnBonds(XR , nb1, 2)
                    J0_up = Invlist_nnBonds(XR , nb1, 3); J0_dn = Invlist_nnBonds(XR , nb1, 4)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case (2)
                    I0_up = Invlist_nnBonds(XR , nb1, 3); I0_dn = Invlist_nnBonds(XR , nb1, 4)
                    J0_up = Invlist_nnBonds(XR , nb1, 1); J0_dn = Invlist_nnBonds(XR , nb1, 2)
                    Zsgn = cmplx( 1.d0,0.d0,Kind(0.d0))
              case default
                 Write(6,*) ' Error in Kin obs '
              end select
              Delta_1 = cmplx(0.d0,0.d0,Kind(0.d0))
              if (I0_up .eq. J0_up) Delta_1 = cmplx(1.d0,0.d0,Kind(0.d0))

              GC_ij(1,1) = Delta_1-GTT(J0_up,I0_up,1); GC_ij(1,2) = -GTT(J0_dn,I0_up,1)
              GC_ij(2,1) = -GTT(J0_up,I0_dn,1); GC_ij(2,2) = Delta_1-GTT(J0_dn,I0_dn,1)

              do ns1 = 1,2
                 Z1 = GC_ij(ns1,ns1)*ZN

                 Z  = Z + Z1*Zsgn
              enddo
           enddo
        endif

        Predefined_Obs_Kin_tau0 = Z
        
     end function Predefined_Obs_Kin_tau0


      
  end submodule ham_QSH_smod
