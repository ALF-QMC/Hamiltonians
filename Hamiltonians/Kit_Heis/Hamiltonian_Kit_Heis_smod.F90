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
!> This File is a template for defining new models. 
!> One can define a new model class by copying this file, replacing alle occurences
!> of ##NAME## by the Hamiltonian name, populating the subroutines below as needed
!> adding the Hamiltonian name to the file Prog/Hamiltonians.list.

!> @details
!> The public variables of this module are the following
!>
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

    submodule (Hamiltonian_main) ham_Kit_Heis_smod

      Use natural_constants, only: pi
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
      
      type, extends(ham_base) :: ham_Kit_Heis
      contains
        ! Set Hamiltonian-specific procedures
        procedure, nopass :: Ham_Set
        procedure, nopass :: Alloc_obs
        procedure, nopass :: Obser
        procedure, nopass :: ObserT
        ! procedure, nopass :: ##PROCEDURE_NAME##  ! Some other procedure defined in ham_base
#ifdef HDF5
        procedure, nopass :: write_parameters_hdf5
#endif
      end type ham_Kit_Heis

      !#PARAMETERS START# VAR_lattice
      Integer            :: L1 = 4   ! Length in direction a_1
      Integer            :: L2 = 4   ! Length in direction a_2
      Character (len=64) :: Lattice_type = 'Honeycomb-Kit'
      Character (len=64) :: Model = 'Kit_Heis'
      !#PARAMETERS END#

      !#PARAMETERS START# VAR_Kit_Heis
      real (Kind=Kind(0.d0)) :: Ham_J      = -0.1d0   ! Nearest neighbor Heisenberg coupling
      real (Kind=Kind(0.d0)) :: Ham_J3     = 0.1d0    ! Third nearest neighbor Heisenberg coupling
      real (Kind=Kind(0.d0)) :: Ham_alphax = -1.d0    ! Kitaev coupling on x-bonds
      real (Kind=Kind(0.d0)) :: Ham_alphay = -1.d0    ! Kitaev coupling on y-bonds
      real (Kind=Kind(0.d0)) :: Ham_alphaz = -1.d0    ! Kitaev coupling on z-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gx     = 0.5d0    ! Gamma coupling on x-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gy     = 0.5d0    ! Gamma coupling on y-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gz     = 0.5d0    ! Gamma coupling on z-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gx_p   = 0.d0     ! Gamma prime coupling on x-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gy_p   = 0.d0     ! Gamma prime coupling on y-bonds
      real (Kind=Kind(0.d0)) :: Ham_Gz_p   = 0.d0     ! Gamma prime coupling on z-bonds
      real (Kind=Kind(0.d0)) :: Hab        = 0.0579d0 ! mu_B * B / K ; Magnetic field strength
      real (Kind=Kind(0.d0)) :: Htheta     = 0.45d0   ! Magnetic field angle (in units of pi)
      real (Kind=Kind(0.d0)) :: Ham_U      = 7.14d0   ! Hubbard U, for freezing out charge fluctuations
      real (Kind=Kind(0.d0)) :: Dtau       = 0.1d0    ! Imaginary time step. Ltrot=Beta/dtau
      real (Kind=Kind(0.d0)) :: Beta       = 5.d0     ! Inverse temperature
      real (Kind=Kind(0.d0)) :: Theta      = 0.d0     ! Projection parameter
      ! logical :: Projector               = .false.  ! Whether the projective algorithm is used
      !#PARAMETERS END#

      real (Kind=Kind(0.d0)) :: Hx, Hy, Hz
      real (Kind=Kind(0.d0)), parameter :: precision = 1.d-12
      real (Kind=Kind(0.d0)), parameter :: ga=2.3d0, gb=2.3d0, gc=1.3d0

      Type (Lattice),       target :: Latt
      Type (Unit_cell),     target :: Latt_unit, latt_unit_2
      Type (Hopping_Matrix_type), Allocatable :: Hopping_Matrix(:)
      Integer, allocatable :: List(:,:), Invlist(:,:)  ! For orbital structure of Unit cell

    contains
      
      module Subroutine Ham_Alloc_Kit_Heis
        allocate(ham_Kit_Heis::ham)
      end Subroutine Ham_Alloc_Kit_Heis

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_Kit_Heis_read_write_parameters.F90"

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
          ! Hard-coded parameters
          N_SUN = 1
          N_FL  = 1
          Symm  = .false.

          ! From dynamically generated file "Hamiltonian_Kit_Heis_read_write_parameters.F90"
          call read_parameters()

          if ( model /= 'Kit_Heis' ) then
             write(error_unit,*) 'Error in Ham_Set: model name does not match the Hamiltonian module'
             CALL Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
          endif
          If (Lattice_type .ne. "Honeycomb-Kit")   then
             Write(error_unit,*) 'Error in Ham_Set: Kit_Heis only implemented on Honeycomb-Kit lattice'
             CALL Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
          Endif

          Ltrot = nint(beta/dtau)
          Thtrot = 0
          if (Projector) Thtrot = nint(theta/dtau)
          Ltrot = Ltrot+2*Thtrot

          ! Setup the Bravais lattice
          call Ham_Latt

          ! Setup single-particle part -> Zeeman term
          call Ham_Hop

          ! Setup the interaction.
          call Ham_V

          ! Setup the trial wave function, in case of a projector approach
          if (Projector) Call Ham_Trial()

#ifdef MPI
          If (Irank_g == 0) then
#endif
             File_info = "info"
#if defined(TEMPERING)
             write(File_info,'(A,I0,A)') "Temp_",igroup,"/info"
#endif
             Open(newunit=unit_info, file=file_info, status="unknown", position="append")
             Write(unit_info,*) '====================================='
             Write(unit_info,*) 'Model is      :  Kit_Heis'
             Write(unit_info,*) 'Lattice is    : ', Lattice_type
             Write(unit_info,*) '# of orbitals : ', Ndim
!              Write(unit_info,*) '# unit cells  : ', Latt%N 
!              Write(unit_info,*) '# of orbitals : ', Latt_unit%Norb
             if (Projector) then
                Write(unit_info,*) 'Projective version'
                Write(unit_info,*) 'Theta         : ', Theta
                Write(unit_info,*) 'Tau_max       : ', beta
                Do nf = 1,N_FL
                   Write(unit_info,*) 'Degen of right trial wave function: ', WF_R(nf)%Degen
                   Write(unit_info,*) 'Degen of left  trial wave function: ', WF_L(nf)%Degen
                enddo
             else
                Write(unit_info,*) 'Finite temperature version'
                Write(unit_info,*) 'Beta          : ', Beta
             endif
             Write(unit_info,*) 'dtau,Ltrot_eff: ', dtau,Ltrot
             Write(unit_info,*) 'N_SUN         : ', N_SUN
             Write(unit_info,*) 'N_FL          : ', N_FL
             Write(unit_info,*) 'Ham_U         : ', Ham_U
             Write(unit_info,*) 'Ham_J         : ', Ham_J
             Write(unit_info,*) 'Ham_J3        : ', Ham_J3
             Write(unit_info,*) 'Ham_alphax    : ', Ham_alphax
             Write(unit_info,*) 'Ham_alphay    : ', Ham_alphay
             Write(unit_info,*) 'Ham_alphaz    : ', Ham_alphaz
             Write(unit_info,*) 'Ham_Gx        : ', Ham_Gx
             Write(unit_info,*) 'Ham_Gy        : ', Ham_Gy
             Write(unit_info,*) 'Ham_Gz        : ', Ham_Gz
             Write(unit_info,*) 'Ham_Gx_p      : ', Ham_Gx_p
             Write(unit_info,*) 'Ham_Gy_p      : ', Ham_Gy_p
             Write(unit_info,*) 'Ham_Gz_p      : ', Ham_Gz_p

             Write(unit_info,*) 'Hab           : ', Hab
             Write(unit_info,*) 'Htheta        : ', Htheta
             Write(unit_info,*) 'Hx            : ', Hx
             Write(unit_info,*) 'Hy            : ', Hy
             Write(unit_info,*) 'Hz            : ', Hz
             Close(unit_info)
#ifdef MPI
          Endif
#endif
        end Subroutine Ham_Set

        subroutine Ham_Latt
          implicit none
          integer :: I, no, nc
          Real (Kind=Kind(0.d0))  :: a1_p(2), a2_p(2), L1_p(2), L2_p(2)

          If (L1==1 .or. L2==1 ) then
            Write(6,*) 'For one-dimensional lattices set : L2 = 1'
            call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
          endif
          Latt_Unit%Norb    = 4
          Latt_Unit%N_coord = 3
          a1_p(1) =  1.D0   ; a1_p(2) =  0.d0
          a2_p(1) =  0.5D0  ; a2_p(2) =  sqrt(3.D0)/2.D0
          Allocate (Latt_Unit%Orb_pos_p(4, 2))
          Latt_Unit%Orb_pos_p(1,:) = 0.d0
          Latt_Unit%Orb_pos_p(2,:) = (a2_p(:) - 0.5D0*a1_p(:) ) * 2.D0/3.D0
          Latt_Unit%Orb_pos_p(3,:) = 0.d0
          Latt_Unit%Orb_pos_p(4,:) = (a2_p(:) - 0.5D0*a1_p(:) ) * 2.D0/3.D0
          L1_p    =  dble(L1) * a1_p
          L2_p    =  dble(L2) * a2_p
          Call Make_Lattice( L1_p, L2_p, a1_p,  a2_p, Latt, nnlist_range_in=2 )
          Ndim = Latt%N*Latt_Unit%Norb
          Allocate (List(Ndim,2), Invlist(Latt%N,Latt_Unit%Norb))
          nc = 0
          Do I = 1,Latt%N
            Do no = 1,Latt_Unit%Norb
                ! For the Honeycomb and pi-flux lattices no = 1,2 corresponds to the A,and B sublattice.
                  ! For  Kit_Heis:      no=1    A, up
                  ! For  Kit_Heis:      no=2    B, up
                  ! For  Kit_Heis:      no=3    A, do
                  ! For  Kit_Heis:      no=4    B, do
                nc = nc + 1
                List(nc,1) = I
                List(nc,2) = no
                Invlist(I,no) = nc 
            Enddo
          Enddo

          Latt_Unit_2%Norb    = 2
          Latt_Unit_2%N_coord = 3
          Allocate (Latt_Unit_2%Orb_pos_p(2,2))
          Latt_Unit_2%Orb_pos_p(1,:) = 0.d0
          Latt_Unit_2%Orb_pos_p(2,:) = (a2_p(:) - 0.5D0*a1_p(:) ) * 2.D0/3.D0
        end subroutine Ham_Latt

        subroutine Ham_hop
          implicit none
          integer :: n, I, no, I1_up, I1_do, I2_up, I2_do

          !ac plane
          Hx=     Hab*sin(Htheta*pi)*ga/(6d0)**0.5d0+Hab*cos(Htheta*pi)*gc/(3d0)**0.5d0
          Hy=     Hab*sin(Htheta*pi)*ga/(6d0)**0.5d0+Hab*cos(Htheta*pi)*gc/(3d0)**0.5d0
          Hz=-2d0*Hab*sin(Htheta*pi)*ga/(6d0)**0.5d0+Hab*cos(Htheta*pi)*gc/(3d0)**0.5d0

          allocate(Op_T(1,N_FL))
          do n = 1,N_FL
            Call Op_make(Op_T(1,n),Ndim)
            DO I = 1, Latt%N
              I1_up = Invlist(I,1) !fAu
              I1_do = Invlist(I,3) !fAd
              I2_up = invlist(I,2) !fBu
              I2_do = invlist(I,4) !fBd
              
              Op_T(1,n)%O(I1_up ,I1_up) = cmplx( 0.5d0*Hz,  0.d0,       kind(0.D0))
              Op_T(1,n)%O(I1_up ,I1_do) = cmplx( 0.5d0*Hx, -0.5d0*Hy,   kind(0.D0))
              Op_T(1,n)%O(I1_up ,I2_up) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I1_up ,I2_do) = cmplx( 0.d0,      0.d0,       kind(0.D0))

              Op_T(1,n)%O(I1_do ,I1_up) = cmplx( 0.5d0*Hx,  0.5d0*Hy,   kind(0.D0))
              Op_T(1,n)%O(I1_do ,I1_do) = cmplx(-0.5d0*Hz,  0.d0,       kind(0.D0))
              Op_T(1,n)%O(I1_do ,I2_up) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I1_do ,I2_do) = cmplx( 0.d0,      0.d0,       kind(0.D0))

              Op_T(1,n)%O(I2_up ,I1_up) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I2_up ,I1_do) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I2_up ,I2_up) = cmplx( 0.5d0*Hz,  0.d0,       kind(0.D0))
              Op_T(1,n)%O(I2_up ,I2_do) = cmplx( 0.5d0*Hx, -0.5d0*Hy,   kind(0.D0))

              Op_T(1,n)%O(I2_do ,I1_up) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I2_do ,I1_do) = cmplx( 0.d0,      0.d0,       kind(0.D0))
              Op_T(1,n)%O(I2_do ,I2_up) = cmplx( 0.5d0*Hx,  0.5d0*Hy,   kind(0.D0))
              Op_T(1,n)%O(I2_do ,I2_do) = cmplx(-0.5d0*Hz,  0.d0,       kind(0.D0))
            Enddo
            Do I = 1,Ndim
                Op_T(1,n)%P(I) = I 
            Enddo
            if ( abs(Hab) < 1.E-6 ) then 
                Op_T(1,n)%g = 0.d0
            else
                Op_T(1,n)%g = -Dtau
            endif
            Op_T(1,n)%alpha=cmplx(0.d0,0.d0, kind(0.D0))
            Call Op_set(Op_T(1,n))
          Enddo
        end subroutine Ham_hop

        subroutine Ham_Trial
          implicit none
            write(error_unit, *) 'Projective algorithm not yet implemented for Kit_Heis model.'
            call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
            !  N_part = Ndim/2
            !  Call Predefined_TrialWaveFunction(Lattice_type ,Ndim,  List,Invlist,Latt, Latt_unit, &
            !       &                            N_part, N_FL,  WF_L, WF_R)
        end subroutine Ham_Trial

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Sets a spin-spin interaction term of the form \f$ J S_a S_b \f$  in the operator format.
!> Here, a and b can be x,y or z.
!> This formulation assumes that all charges fluctuations are frozen out, e.g. by a large Hubbard U.
!> Thereby, \f$ c^\dagger_{i,\uparrow} c_{i,\uparrow} + c^\dagger_{i,\downarrow} c_{i,\downarrow} = 1 \f$.
!> The spin operators are represented in terms of fermionic operators as follows:
!> \f$ S^{\alpha}_i = \frac{1}{2}(c^\dagger_i \sigma^\alpha c_i)\f$
!> Here, \f$ c^\dagger_i = (c^\dagger_{i,\uparrow}, c^\dagger_{i,\downarrow}) \f$ is a spinor of fermionic creation operators and \f$ \sigma^\alpha \f$ are the Pauli matrices.
!> The exact term implemented here is \f$g \abs{J} / 8 (S^a_i + conjg(g) J/\abs{J} S^b_j)^2 = S^a_i S^b_j + const.\f$.
!> The term includes a gauge factor g = gauge, which can possibly be used to alleviate the sign problem.
!> The gauge factor is a pure phase factor, defined as \f$ g = e^{i \phi} \f$ with \f$ \phi = \text{gauge\_phi} \pi \f$.
!--------------------------------------------------------------------
        subroutine Predefined_Int_SaSb(Op, a, b, I1_up, I2_up, I1_do, I2_do, J, dtau, gauge_phi)
          implicit none
          type(Operator), intent(out) :: Op
          Character (len=1), intent(in) :: a, b
          integer, intent(in) :: I1_up, I2_up, I1_do, I2_do
          real(Kind=Kind(0.d0)), intent(in) :: J, dtau
          real(Kind=Kind(0.d0)), intent(in) :: gauge_phi

          real(Kind=Kind(0.d0)), parameter :: pi = acos(-1.d0)
          real(Kind=Kind(0.d0)) :: J_sgn
          complex(Kind=Kind(0.d0)) :: gauge

            gauge = cmplx(cos(gauge_phi*pi), sin(gauge_phi*pi), kind(0.D0))
            if (abs(abs(gauge)-1.d0) > 1.d-8) then
                write(error_unit, *) 'Error in Predefined_Int_SaSb: gauge must be a pure phase factor with absolute value 1.'
                call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
            endif
            if (abs(J).le.1D-12)then
              J_sgn = 0.d0
            else
              J_sgn = J/abs(J)
            end if

            call Op_make(Op,4)
            Op%P(1) = I1_up
            Op%P(2) = I2_up
            Op%P(3) = I1_do
            Op%P(4) = I2_do
            select case (a)
            case('x')
               Op%O(1,3) = cmplx(1.d0, 0.d0, kind(0.D0))
               Op%O(3,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
            case('y')
               Op%O(1,3) = cmplx(0.d0, -1.d0, kind(0.D0))
               Op%O(3,1) = cmplx(0.d0 ,1.d0, kind(0.D0))
            case('z')
               Op%O(1,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
               Op%O(3,3) = cmplx(-1.d0 ,0.d0, kind(0.D0))
            end select
            select case (b)
            case('x')
               Op%O(2,4) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
               Op%O(4,2) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
            case('y')
               Op%O(2,4) = conjg(gauge)*cmplx(0.d0, -J_sgn*1.d0, kind(0.D0))
               Op%O(4,2) = conjg(gauge)*cmplx(0.d0 ,J_sgn*1.d0, kind(0.D0))
            case('z')
               Op%O(2,2) = conjg(gauge)*cmplx(J_sgn*1.d0 ,0.d0, kind(0.D0))
               Op%O(4,4) = conjg(gauge)*cmplx(-J_sgn*1.d0 ,0.d0, kind(0.D0))
            end select
            Op%g     = SQRT(gauge*CMPLX(-dtau*abs(J)/8d0, 0.D0, kind(0.D0)))
            Op%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
            Op%type  = 2
            call Op_set(Op)
        end subroutine Predefined_Int_SaSb

        subroutine add_int_SaSb(Ops, nc, a, b, I1_up, I2_up, I1_do, I2_do, J, dtau, gauge_phi)
          implicit none
          type(Operator), intent(inout) :: Ops(:, :)
          integer, intent(inout) :: nc
          Character (len=1), intent(in) :: a, b
          integer, intent(in) :: I1_up, I2_up, I1_do, I2_do
          real(Kind=Kind(0.d0)), intent(in) :: J, dtau
          real(Kind=Kind(0.d0)), intent(in) :: gauge_phi

          if (abs(J) > precision) then
            nc = nc + 1
            call Predefined_Int_SaSb(Ops(nc, 1), a, b, I1_up, I2_up, I1_do, I2_do, J, dtau, gauge_phi)
          endif
        end subroutine add_int_SaSb

!--------------------------------------------------------------------
!> @author 
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
        Subroutine Ham_V
          
          Implicit none 
          
          Integer :: I, nc, nc1
          Integer :: I1_up, I1_do, I2_up, I2_do
          Integer :: N_ops

          N_ops = 0
          if (abs(Ham_J) > precision) N_ops = N_ops + 2*Latt%N*Latt_unit%N_coord   ! Heisenberg term, real and imaginary part
          if (abs(Ham_J3) > precision) N_ops = N_ops + 3*Latt%N*Latt_unit%N_coord  ! Third nearest neighbor Heisenberg term
          if (abs(Ham_alphax) > precision) N_ops = N_ops + Latt%N                  ! Kitaev term, x-bonds
          if (abs(Ham_alphay) > precision) N_ops = N_ops + Latt%N                  ! Kitaev term, y-bonds
          if (abs(Ham_alphaz) > precision) N_ops = N_ops + Latt%N                  ! Kitaev term, z-bonds
          if (abs(Ham_Gx) > precision) N_ops = N_ops + 2*Latt%N                    ! Gamma term, x-bonds
          if (abs(Ham_Gy) > precision) N_ops = N_ops + 2*Latt%N                    ! Gamma term, y-bonds
          if (abs(Ham_Gz) > precision) N_ops = N_ops + 2*Latt%N                    ! Gamma term, z-bonds
          if (abs(Ham_Gx_p) > precision) N_ops = N_ops + 4*Latt%N                  ! Gamma prime term, x-bonds
          if (abs(Ham_Gy_p) > precision) N_ops = N_ops + 4*Latt%N                  ! Gamma prime term, y-bonds
          if (abs(Ham_Gz_p) > precision) N_ops = N_ops + 4*Latt%N                  ! Gamma prime term, z-bonds
          N_ops = N_ops + 2*Latt%N                                                 ! Hubbard U term
          Allocate(Op_V(N_ops, N_FL))

          nc = 0
          if (abs(Ham_J) > precision) then
          do I = 1,Latt%N  ! Heisenberg, real part
            do nc1 = 1, Latt_unit%N_coord
                I1_up = Invlist(I,1)
                I1_do = Invlist(I,3)
                select case(nc1)
                case(1)
                  I2_up = invlist(I,2)
                  I2_do = invlist(I,4)
                case(2)
                  I2_up = invlist(Latt%nnlist(I,1, -1),2)
                  I2_do = invlist(Latt%nnlist(I,1, -1),4)
                case(3)
                  I2_up = invlist(Latt%nnlist(I,0, -1),2)
                  I2_do = invlist(Latt%nnlist(I,0, -1),4)
                end select
                nc = nc + 1
                call Op_make(Op_V(nc,1),4)
                Op_V(nc,1)%P(1) = I1_up
                Op_V(nc,1)%P(2) = I2_up
                Op_V(nc,1)%P(3) = I1_do
                Op_V(nc,1)%P(4) = I2_do
                Op_V(nc,1)%O(1,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(2,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(3,4) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%O(4,3) = cmplx(1.d0 ,0.d0, kind(0.D0))
                Op_V(nc,1)%g     = SQRT(CMPLX(DTAU*Ham_J/8.d0, 0.D0, kind(0.D0)))
                Op_V(nc,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
                Op_V(nc,1)%type  = 2
                Call Op_set( Op_V(nc,1) )
            Enddo
          Enddo
          endif

          if (abs(Ham_J) > precision) then
          do I = 1,Latt%N  ! Heisenberg, imaginary part
            do nc1 = 1, Latt_unit%N_coord
              I1_up = Invlist(I,1)
              I1_do = Invlist(I,3)
              select case(nc1)
              case(1)
                I2_up = invlist(I,2)
                I2_do = invlist(I,4)
              case(2)
                I2_up = invlist(Latt%nnlist(I,1, -1),2)
                I2_do = invlist(Latt%nnlist(I,1, -1),4)
              case(3)
                I2_up = invlist(Latt%nnlist(I,0, -1),2)
                I2_do = invlist(Latt%nnlist(I,0, -1),4)
              end select
              nc = nc + 1
              call Op_make(Op_V(nc,1),4)
              Op_V(nc,1)%P(1) = I1_up
              Op_V(nc,1)%P(2) = I2_up
              Op_V(nc,1)%P(3) = I1_do
              Op_V(nc,1)%P(4) = I2_do
              Op_V(nc,1)%O(1,2) = cmplx(0.d0 , 1.d0, kind(0.D0))
              Op_V(nc,1)%O(2,1) = cmplx(0.d0 ,-1.d0, kind(0.D0))
              Op_V(nc,1)%O(3,4) = cmplx(0.d0 , 1.d0, kind(0.D0))
              Op_V(nc,1)%O(4,3) = cmplx(0.d0 ,-1.d0, kind(0.D0))
              Op_V(nc,1)%g     = SQRT(CMPLX(DTAU*Ham_J/8.d0, 0.D0, kind(0.D0)))
              Op_V(nc,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,1)%type  = 2
              Call Op_set( Op_V(nc,1) )
            Enddo
          Enddo
          endif

          do I = 1,Latt%N  ! J3-SxSx
            do nc1 = 1, Latt_unit%N_coord
              I1_up = Invlist(I,1)
              I1_do = Invlist(I,3)
              select case(nc1)
              case(1)
                  I2_up = invlist(Latt%nnlist(I, 1, 0),2)
                  I2_do = invlist(Latt%nnlist(I, 1, 0),4)
              case(2)
                  I2_up = invlist(Latt%nnlist(I,-1, 0),2)
                  I2_do = invlist(Latt%nnlist(I,-1, 0),4)
              case(3)
                  I2_up = invlist(Latt%nnlist(I, 1, -2),2)
                  I2_do = invlist(Latt%nnlist(I, 1, -2),4)
              end select
              call add_int_SaSb(Op_V, nc,'x', 'x', I1_up, I2_up, I1_do, I2_do, Ham_J3, dtau, 1.d0)
            Enddo
          Enddo

          do I = 1,Latt%N  ! J3-SySy
            do nc1 = 1, Latt_unit%N_coord
              I1_up = Invlist(I,1)
              I1_do = Invlist(I,3)
              select case(nc1)
              case(1)
                  I2_up = invlist(Latt%nnlist(I, 1, 0),2)
                  I2_do = invlist(Latt%nnlist(I, 1, 0),4)
              case(2)
                  I2_up = invlist(Latt%nnlist(I,-1, 0),2)
                  I2_do = invlist(Latt%nnlist(I,-1, 0),4)
              case(3)
                  I2_up = invlist(Latt%nnlist(I, 1, -2),2)
                  I2_do = invlist(Latt%nnlist(I, 1, -2),4)
              end select
              call add_int_SaSb(Op_V, nc,'y', 'y', I1_up, I2_up, I1_do, I2_do, Ham_J3, dtau, 1.d0)
            Enddo
          Enddo

          do I = 1,Latt%N  ! J3-SzSz
            do nc1 = 1, Latt_unit%N_coord
              I1_up = Invlist(I,1)
              I1_do = Invlist(I,3)
              select case(nc1)
              case(1)
                  I2_up = invlist(Latt%nnlist(I, 1, 0),2)
                  I2_do = invlist(Latt%nnlist(I, 1, 0),4)
              case(2)
                  I2_up = invlist(Latt%nnlist(I,-1, 0),2)
                  I2_do = invlist(Latt%nnlist(I,-1, 0),4)
              case(3)
                  I2_up = invlist(Latt%nnlist(I, 1, -2),2)
                  I2_do = invlist(Latt%nnlist(I, 1, -2),4)
              end select
              call add_int_SaSb(Op_V, nc,'z', 'z', I1_up, I2_up, I1_do, I2_do, Ham_J3, dtau, 1.d0)
            Enddo
          Enddo

          do I = 1,Latt%N  ! Kitaev
            do nc1 = 1, Latt_unit%N_coord
              I1_up = Invlist(I,1)
              I1_do = Invlist(I,3)
              select case(nc1)
              case(1)  ! Z-bonds
                I2_up = invlist(I,2)
                I2_do = invlist(I,4)
                ! Kitaev
                call add_int_SaSb(Op_V, nc,'z', 'z', I1_up, I2_up, I1_do, I2_do, Ham_alphaz, dtau, 1.d0)
                ! Gamma terms
                call add_int_SaSb(Op_V, nc,'x', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gz, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'y', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gz, dtau, 1.d0)
                ! Gamma prime terms
                call add_int_SaSb(Op_V, nc,'z', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gz_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'z', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gz_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'y', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gz_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'x', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gz_p, dtau, 1.d0)
              case(2)  ! Y-bonds
                I2_up = invlist(Latt%nnlist(I,1, -1),2)
                I2_do = invlist(Latt%nnlist(I,1, -1),4)
                ! Kitaev
                call add_int_SaSb(Op_V, nc,'y', 'y', I1_up, I2_up, I1_do, I2_do, Ham_alphay, dtau, 1.d0)
                ! Gamma terms
                call add_int_SaSb(Op_V, nc,'z', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gy, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'x', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gy, dtau, 1.d0)
                ! Gamma prime terms
                call add_int_SaSb(Op_V, nc,'y', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gy_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'y', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gy_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'x', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gy_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'z', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gy_p, dtau, 1.d0)
              case(3)  ! X-bonds
                I2_up = invlist(Latt%nnlist(I,0, -1),2)
                I2_do = invlist(Latt%nnlist(I,0, -1),4)
                ! Kitaev
                call add_int_SaSb(Op_V, nc,'x', 'x', I1_up, I2_up, I1_do, I2_do, Ham_alphax, dtau, 1.d0)
                ! Gamma terms
                call add_int_SaSb(Op_V, nc,'y', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gx, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'z', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gx, dtau, 1.d0)
                ! Gamma prime terms
                call add_int_SaSb(Op_V, nc,'x', 'y', I1_up, I2_up, I1_do, I2_do, Ham_Gx_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'x', 'z', I1_up, I2_up, I1_do, I2_do, Ham_Gx_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'z', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gx_p, dtau, 1.d0)
                call add_int_SaSb(Op_V, nc,'y', 'x', I1_up, I2_up, I1_do, I2_do, Ham_Gx_p, dtau, 1.d0)
              end select
            Enddo
          Enddo


          do I = 1,Latt%N  ! Hubbard U term
            do nc1 = 1, 2
              select case(nc1)
              case(1) ! Orbital A
                I1_up = invlist(I,1)
                I1_do = invlist(I,3)
              case(2) ! Orbital B
                I1_up = invlist(I,2)
                I1_do = invlist(I,4)
              end select
              nc = nc + 1
              call Op_make(Op_V(nc,1),2)
              Op_V(nc,1)%P(1) = I1_up
              Op_V(nc,1)%P(2) = I1_do
              ! This is for SU(2) Hubbard
              Op_V(nc,1)%O(1,1) = cmplx(1.d0 ,0.d0, kind(0.D0))
              Op_V(nc,1)%O(2,2) = cmplx(1.d0 ,0.d0, kind(0.D0))
              Op_V(nc,1)%g     = SQRT(CMPLX(-DTAU*Ham_U/2.d0, 0.D0, kind(0.D0)))
              Op_V(nc,1)%alpha = cmplx(-1.d0, 0.d0, kind(0.D0))
              ! This is for M_z Hubbard (may yield a better sign problem)
              ! Op_V(nc,1)%O(1,1) = cmplx( 1.d0 ,0.d0, kind(0.D0))
              ! Op_V(nc,1)%O(2,2) = cmplx(-1.d0 ,0.d0, kind(0.D0))
              ! Op_V(nc,1)%g     = SQRT(CMPLX(DTAU*Ham_U/2.d0, 0.D0, kind(0.D0)))
              ! Op_V(nc,1)%alpha = cmplx(0.d0, 0.d0, kind(0.D0))
              Op_V(nc,1)%type  = 2
              Call Op_set( Op_V(nc,1) )
            Enddo
          Enddo

          if (nc .ne. size(Op_V,1)) then
              write(error_unit, *) 'Error in Ham_V: number of interaction terms doe not fit size of Op_V.'
              write(error_unit, *) 'Number of interaction terms counted: ', nc
              write(error_unit, *) 'Size of Op_V: ', size(Op_V,1)
              call Terminate_on_error(ERROR_HAMILTONIAN, __FILE__, __LINE__)
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
          Integer    ::  I, N, Nt
          Character (len=64) ::  Filename
          Character (len=:), allocatable ::  Channel

          ! Scalar observables
          Allocate ( Obs_scal(24) )
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
              N = 1;   Filename ="Ek"
             case (6)
              N = 1;   Filename ="Eg"
             case (7)
              N = 1;   Filename ="Ej1"
             case (8)
              N = 1;   Filename ="Ej3"
             case (9)
              N = 1;   Filename ="Eh"

           case (10)
              N = 1;   Filename ="EKxx"
           case (11)
              N = 1;   Filename ="EKyy"
           case (12)
              N = 1;   Filename ="EKzz"
           case (13)
              N = 1;   Filename ="EKyz"
           case (14)
              N = 1;   Filename ="EKzx"
           case (15)
              N = 1;   Filename ="EKxy"
           case (16)
              N = 1;   Filename ="EJ1xx"
           case (17)
              N = 1;   Filename ="EJ1yy"
           case (18)
              N = 1;   Filename ="EJ1zz"
           case (19)
              N = 1;   Filename ="EJ3xx"
           case (20)
              N = 1;   Filename ="EJ3yy"
           case (21)
              N = 1;   Filename ="EJ3zz"
           case (22)
              N = 1;   Filename ="Esx"
           case (23)
              N = 1;   Filename ="Esy"
           case (24)
              N = 1;   Filename ="Esz"              
           case default
              Write(6,*) ' Error in Alloc_obs '
             end select
             Call Obser_Vec_make(Obs_scal(I),N,Filename)
           enddo

           Allocate ( Obs_eq(11) )
           Nt = 1
           Channel = '--'
           Do I = 1,Size(Obs_eq,1)
             select case (I)
             case (1)
                Filename ="Green"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             case (2)
                Filename ="Spinxx"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (3)
                Filename ="Spinyy"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (4)
                Filename ="Spinzz"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (5)
                Filename ="Spinxy"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (6)
                Filename ="Spinxz"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (7)
                Filename ="Spinyx"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (8)
                Filename ="Spinyz"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (9)
                Filename ="Spinzx"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (10)
                Filename ="Spinzy"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (11)
                Filename ="Den"
                call Obser_Latt_make(Obs_eq(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
           enddo

           If (Ltau == 1) then
           Allocate ( Obs_tau(11) )
           Nt = Ltrot+1-2*Thtrot
           Do I = 1,Size(Obs_tau,1)
             select case (I)
             case (1)
                Filename ="Green"; Channel = 'P'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
             case (2)
                Filename ="Spinxx"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (3)
                Filename ="Spinyy"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (4)
                Filename ="Spinzz"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (5)
                Filename ="Spinxy"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (6)
                Filename ="Spinxz"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (7)
                Filename ="Spinyx"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (8)
                Filename ="Spinyz"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (9)
                Filename ="Spinzx"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (10)
                Filename ="Spinzy"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case (11)
                Filename ="Den"; Channel = 'PH'
                call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit_2, Channel, dtau)
             case default
                Write(6,*) ' Error in Alloc_obs '
             end select
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
        subroutine Obser(GR,Phase,Ntau, Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
          Integer,                   INTENT(IN) :: Ntau
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

          !Local
          Complex (Kind=Kind(0.d0)), allocatable :: GRC(:,:,:)
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          Integer :: I, J, nf
          ! Add local variables as needed
          Complex (Kind=Kind(0.d0)) :: Zrho, Zkin, ZPot, Z
          Integer :: imj, I1, J1, no_I, no_J,n

          Integer :: I1_up, I1_do, I2_up, I2_do, nc1, Norb
          Complex (Kind=Kind(0.d0)) :: EKxx,EKyy,EKzz,EKyz,EKzx,EKxy
          Complex (Kind=Kind(0.d0)) :: EJ1xx,EJ1yy,EJ1zz,EJ3xx,EJ3yy,EJ3zz
          Complex (Kind=Kind(0.d0)) :: Esx,Esy,Esz

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

          ZS = ZS*Mc_step_weight

          allocate(GRC(Ndim,Ndim,N_FL))
          
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
          Do nf = 1,N_FL
             Do I = 1,Latt%N
                Do n = 1,Latt_unit%N_coord
                   I1 = invlist(I,1)
                   If (n == 1)  J1 = invlist(I,2)
                   If (n == 2)  J1 = invlist(Latt%nnlist(I,1,-1),2)
                   If (n == 3)  J1 = invlist(Latt%nnlist(I,0,-1),2)
                   Zkin = Zkin +  Grc( I1,J1, nf ) + Grc(J1,I1,nf)
                Enddo
             Enddo
          Enddo
          Zkin = -Zkin * dble(N_SUN)
          Obs_scal(1)%Obs_vec(1)  =    Obs_scal(1)%Obs_vec(1) + Zkin *ZP* ZS


          ZPot = cmplx(0.d0, 0.d0, kind(0.D0))
             do I = 1,Latt%N
                do nc1 = 1, 2
                   select case(nc1)
                   case(1) ! Orbital A
                      I1_up = invlist(I,1)
                      I1_do = invlist(I,3)
                   case(2) ! Orbital B 
                      I1_up = invlist(I,2)
                      I1_do = invlist(I,4)
                   end select
                   ZPot = ZPot + Grc(I1_up,I1_up,1) * Grc(I1_do,I1_do,1)+Grc(I1_up,I1_do,1)* Gr(I1_up,I1_do,1)
                Enddo
             Enddo
            Zpot = Zpot/dble(Latt%N*2)

          Obs_scal(2)%Obs_vec(1)  =  Obs_scal(2)%Obs_vec(1) + Zpot * ZP*ZS


          Zrho = cmplx(0.d0,0.d0, kind(0.D0))
          Do nf = 1,N_FL
             Do I = 1,Ndim
                Zrho = Zrho + Grc(i,i,nf) 
             enddo
          enddo
          Zrho = Zrho* dble(N_SUN)
          Obs_scal(3)%Obs_vec(1)  =    Obs_scal(3)%Obs_vec(1) + Zrho * ZP*ZS
          Obs_scal(4)%Obs_vec(1)  =    Obs_scal(4)%Obs_vec(1) + (Zkin + Zpot)*ZP*ZS


          ! Compute equal-time correlations
          ! Compute spin-spin, Green, and den-den correlation functions  !  This is general N_SUN, and  N_FL = 1
          DO I = 1,Size(Obs_eq,1)
             Obs_eq(I)%N        = Obs_eq(I)%N + 1
             Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
          ENDDO

            Norb = Latt_unit%Norb
            Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))

            EKxx=cmplx(0.d0, 0.d0, kind(0.D0))
            EKyy=cmplx(0.d0, 0.d0, kind(0.D0))
            EKzz=cmplx(0.d0, 0.d0, kind(0.D0))
            EKyz=cmplx(0.d0, 0.d0, kind(0.D0))
            EKzx=cmplx(0.d0, 0.d0, kind(0.D0))
            EKxy=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ1xx=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ1yy=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ1zz=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ3xx=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ3yy=cmplx(0.d0, 0.d0, kind(0.D0))
            EJ3zz=cmplx(0.d0, 0.d0, kind(0.D0))
            Esx=cmplx(0.d0, 0.d0, kind(0.D0))
            Esy=cmplx(0.d0, 0.d0, kind(0.D0))
            Esz=cmplx(0.d0, 0.d0, kind(0.D0))
            

            do I = 1,Latt%N  !  Kitaev
               do nc1 = 1, Latt_unit%N_coord
                  I1_up = Invlist(I,1)
                  I1_do = Invlist(I,3)
                  select case(nc1)

                  case(1) ! X-bonds X-bonds-yz                     
                     I2_up = invlist(Latt%nnlist(I,0, -1),2)
                     I2_do = invlist(Latt%nnlist(I,0, -1),4)
                     
                     !xx 4*
                     EKxx=Ekxx+( &
                          & GRC(I1_up, I1_do,1)*GRC(I2_up, I2_do,1)+GRC(I1_up, I2_do,1)*GR(I1_do, I2_up,1) &
                          &+GRC(I1_up, I1_do,1)*GRC(I2_do, I2_up,1)+GRC(I1_up, I2_up,1)*GR(I1_do, I2_do,1) &
                          &+GRC(I1_do, I1_up,1)*GRC(I2_up, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_up, I2_up,1) &
                          &+GRC(I1_do, I1_up,1)*GRC(I2_do, I2_up,1)+GRC(I1_do, I2_up,1)*GR(I1_up, I2_do,1) )

                     !yz i4*
                     EKyz=EKyz+( &
                          & GRC(I1_up, I1_do,1)*GRC(I2_up , I2_up ,1)+GRC(I1_up , I2_up ,1)*GR(I1_do, I2_up ,1) &
                          &-GRC(I1_up, I1_do,1)*GRC(I2_do, I2_do,1)-GRC(I1_up , I2_do,1)*GR(I1_do, I2_do,1) &
                          &-GRC(I1_do, I1_up,1)*GRC(I2_up , I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_up , I2_up ,1) &
                          &+GRC(I1_do, I1_up,1)*GRC(I2_do, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_up , I2_do,1) )

                  case(2) ! Y-bonds Y-bonds-zx
                     I2_up = invlist(Latt%nnlist(I,1, -1),2)
                     I2_do = invlist(Latt%nnlist(I,1, -1),4)

                     !yy (-4)*
                     EKyy=Ekyy+( &
                          & GRC(I1_up , I1_do,1)*GRC(I2_up , I2_do,1)+GRC(I1_up , I2_do,1)*GR(I1_do, I2_up ,1) &
                          &-GRC(I1_up , I1_do,1)*GRC(I2_do, I2_up ,1)-GRC(I1_up , I2_up ,1)*GR(I1_do, I2_do,1) &
                          &-GRC(I1_do, I1_up ,1)*GRC(I2_up , I2_do,1)-GRC(I1_do, I2_do,1)*GR(I1_up , I2_up ,1) &
                          &+GRC(I1_do, I1_up ,1)*GRC(I2_do, I2_up ,1)+GRC(I1_do, I2_up ,1)*GR(I1_up , I2_do,1) )

                     !zx 4*
                     EKzx=EKzx+( &
                          & GRC(I1_up , I1_up ,1)*GRC(I2_up , I2_do,1)+GRC(I1_up , I2_do,1)*GR(I1_up , I2_up ,1) &
                          &+GRC(I1_up , I1_up ,1)*GRC(I2_do, I2_up ,1)+GRC(I1_up , I2_up ,1)*GR(I1_up , I2_do,1) &
                          &-GRC(I1_do, I1_do,1)*GRC(I2_up , I2_do,1)-GRC(I1_do, I2_do,1)*GR(I1_do, I2_up ,1) &
                          &-GRC(I1_do, I1_do,1)*GRC(I2_do, I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_do, I2_do,1) )
                     
                  case(3) ! Z-bonds  Z-bonds-xy
                     I2_up = invlist(I,2)
                     I2_do = invlist(I,4)
                     
                     !zz 4*
                     Ekzz=Ekzz+( &
                          & GRC(I1_up , I1_up ,1)*GRC(I2_up , I2_up ,1)+GRC(I1_up , I2_up ,1)*GR(I1_up , I2_up ,1) &
                          &-GRC(I1_up , I1_up ,1)*GRC(I2_do, I2_do,1)-GRC(I1_up , I2_do,1)*GR(I1_up , I2_do,1) &
                          &-GRC(I1_do, I1_do,1)*GRC(I2_up , I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_do, I2_up ,1) &
                          &+GRC(I1_do, I1_do,1)*GRC(I2_do, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_do, I2_do,1) )

                     !xy i4*
                     EKxy=EKxy+( &
                          & GRC(I1_up , I1_do,1)*GRC(I2_up , I2_do,1)+GRC(I1_up , I2_do,1)*GR(I1_do, I2_up ,1) &
                          &-GRC(I1_up , I1_do,1)*GRC(I2_do, I2_up ,1)-GRC(I1_up , I2_up ,1)*GR(I1_do, I2_do,1) &
                          &+GRC(I1_do, I1_up ,1)*GRC(I2_up , I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_up , I2_up ,1) &
                          &-GRC(I1_do, I1_up ,1)*GRC(I2_do, I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_up , I2_do,1) )

                  end select
               Enddo
            Enddo
            Obs_scal(5)%Obs_vec(1)  =  Obs_scal(5)%Obs_vec(1) + dble(Ham_alphax)/4d0*(EKxx - Ekyy + Ekzz) * ZP*ZS
            Obs_scal(6)%Obs_vec(1)  =  Obs_scal(6)%Obs_vec(1) + dble(Ham_Gx)/4d0*ZP*ZS*(&
                                        cmplx(aimag(EKyz),-real(EKyz),kind(0.D0))&
                                       +EKzx&
                                       +cmplx(aimag(EKxy),-real(EKxy),kind(0.D0))&
                                        )
            Obs_scal(10)%Obs_vec(1)  =  Obs_scal(10)%Obs_vec(1) + 1d0/4d0*EKxx * ZP*ZS
            Obs_scal(11)%Obs_vec(1)  =  Obs_scal(11)%Obs_vec(1) + 1d0/4d0*EKyy * ZP*ZS
            Obs_scal(12)%Obs_vec(1)  =  Obs_scal(12)%Obs_vec(1) + 1d0/4d0*EKzz * ZP*ZS
            Obs_scal(13)%Obs_vec(1)  =  Obs_scal(13)%Obs_vec(1) + 1d0/4d0*EKyz * ZP*ZS
            Obs_scal(14)%Obs_vec(1)  =  Obs_scal(14)%Obs_vec(1) + 1d0/4d0*EKzx * ZP*ZS
            Obs_scal(15)%Obs_vec(1)  =  Obs_scal(15)%Obs_vec(1) + 1d0/4d0*EKxy * ZP*ZS
            
            
            do I = 1,Latt%N !   Heisenberg
               do nc1 = 1, Latt_unit%N_coord
                  I1_up = Invlist(I,1)
                  I1_do = Invlist(I,3)
                  select case(nc1)
                  case(1)
                     I2_up = invlist(I,2)
                     I2_do = invlist(I,4)
                  case(2)
                     I2_up = invlist(Latt%nnlist(I,1, -1),2)
                     I2_do = invlist(Latt%nnlist(I,1, -1),4)
                  case(3)
                     I2_up = invlist(Latt%nnlist(I,0, -1),2)
                     I2_do = invlist(Latt%nnlist(I,0, -1),4)
                  end select

                  EJ1xx=EJ1xx+( &
                       & GRC(I1_up, I1_do,1)*GRC(I2_up, I2_do,1)+GRC(I1_up, I2_do,1)*GR(I1_do, I2_up,1) &
                       &+GRC(I1_up, I1_do,1)*GRC(I2_do, I2_up,1)+GRC(I1_up, I2_up,1)*GR(I1_do, I2_do,1) &
                       &+GRC(I1_do, I1_up,1)*GRC(I2_up, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_up, I2_up,1) &
                       &+GRC(I1_do, I1_up,1)*GRC(I2_do, I2_up,1)+GRC(I1_do, I2_up,1)*GR(I1_up, I2_do,1) )

                  EJ1yy=EJ1yy+( &
                       & GRC(I1_up , I1_do,1)*GRC(I2_up , I2_do,1)+GRC(I1_up , I2_do,1)*GR(I1_do, I2_up ,1) &
                       &-GRC(I1_up , I1_do,1)*GRC(I2_do, I2_up ,1)-GRC(I1_up , I2_up ,1)*GR(I1_do, I2_do,1) &
                       &-GRC(I1_do, I1_up ,1)*GRC(I2_up , I2_do,1)-GRC(I1_do, I2_do,1)*GR(I1_up , I2_up ,1) &
                       &+GRC(I1_do, I1_up ,1)*GRC(I2_do, I2_up ,1)+GRC(I1_do, I2_up ,1)*GR(I1_up , I2_do,1) )

                  EJ1zz=EJ1zz+( &
                       & GRC(I1_up , I1_up ,1)*GRC(I2_up , I2_up ,1)+GRC(I1_up , I2_up ,1)*GR(I1_up , I2_up ,1) &
                       &-GRC(I1_up , I1_up ,1)*GRC(I2_do, I2_do,1)-GRC(I1_up , I2_do,1)*GR(I1_up , I2_do,1) &
                       &-GRC(I1_do, I1_do,1)*GRC(I2_up , I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_do, I2_up ,1) &
                       &+GRC(I1_do, I1_do,1)*GRC(I2_do, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_do, I2_do,1) )
                  
               Enddo
            Enddo
            Obs_scal(7)%Obs_vec(1)  =  Obs_scal(7)%Obs_vec(1) + dble(Ham_J)/4d0*(EJ1xx - EJ1yy + EJ1zz) * ZP*ZS

            Obs_scal(16)%Obs_vec(1)  =  Obs_scal(16)%Obs_vec(1) + 1d0/4d0*EJ1xx * ZP*ZS
            Obs_scal(17)%Obs_vec(1)  =  Obs_scal(17)%Obs_vec(1) + 1d0/4d0*EJ1yy * ZP*ZS
            Obs_scal(18)%Obs_vec(1)  =  Obs_scal(18)%Obs_vec(1) + 1d0/4d0*EJ1zz * ZP*ZS

            
            do I = 1,Latt%N !  J3-SxSx
               do nc1 = 1, Latt_unit%N_coord
                  I1_up = Invlist(I,1)
                  I1_do = Invlist(I,3)
                  select case(nc1)
                  case(1)
                     I2_up = invlist(Latt%nnlist(I, 1, 0),2)
                     I2_do = invlist(Latt%nnlist(I, 1, 0),4)
                  case(2)
                     I2_up = invlist(Latt%nnlist(I,-1, 0),2)
                     I2_do = invlist(Latt%nnlist(I,-1, 0),4)
                  case(3)
                     I2_up = invlist(Latt%nnlist(I, 1, -2),2)
                     I2_do = invlist(Latt%nnlist(I, 1, -2),4)
                  end select

                  EJ3xx=EJ3xx+( &
                       & GRC(I1_up, I1_do,1)*GRC(I2_up, I2_do,1)+GRC(I1_up, I2_do,1)*GR(I1_do, I2_up,1) &
                       &+GRC(I1_up, I1_do,1)*GRC(I2_do, I2_up,1)+GRC(I1_up, I2_up,1)*GR(I1_do, I2_do,1) &
                       &+GRC(I1_do, I1_up,1)*GRC(I2_up, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_up, I2_up,1) &
                       &+GRC(I1_do, I1_up,1)*GRC(I2_do, I2_up,1)+GRC(I1_do, I2_up,1)*GR(I1_up, I2_do,1) )

                  EJ3yy=EJ3yy+( &
                       & GRC(I1_up , I1_do,1)*GRC(I2_up , I2_do,1)+GRC(I1_up , I2_do,1)*GR(I1_do, I2_up ,1) &
                       &-GRC(I1_up , I1_do,1)*GRC(I2_do, I2_up ,1)-GRC(I1_up , I2_up ,1)*GR(I1_do, I2_do,1) &
                       &-GRC(I1_do, I1_up ,1)*GRC(I2_up , I2_do,1)-GRC(I1_do, I2_do,1)*GR(I1_up , I2_up ,1) &
                       &+GRC(I1_do, I1_up ,1)*GRC(I2_do, I2_up ,1)+GRC(I1_do, I2_up ,1)*GR(I1_up , I2_do,1) )

                  EJ3zz=EJ3zz+( &
                       & GRC(I1_up , I1_up ,1)*GRC(I2_up , I2_up ,1)+GRC(I1_up , I2_up ,1)*GR(I1_up , I2_up ,1) &
                       &-GRC(I1_up , I1_up ,1)*GRC(I2_do, I2_do,1)-GRC(I1_up , I2_do,1)*GR(I1_up , I2_do,1) &
                       &-GRC(I1_do, I1_do,1)*GRC(I2_up , I2_up ,1)-GRC(I1_do, I2_up ,1)*GR(I1_do, I2_up ,1) &
                       &+GRC(I1_do, I1_do,1)*GRC(I2_do, I2_do,1)+GRC(I1_do, I2_do,1)*GR(I1_do, I2_do,1) )                  

               Enddo
            Enddo
            Obs_scal(8)%Obs_vec(1)  =  Obs_scal(8)%Obs_vec(1) + dble(Ham_J3)/4d0*(EJ3xx - EJ3yy + EJ3zz) * ZP*ZS
            
            Obs_scal(19)%Obs_vec(1)  =  Obs_scal(19)%Obs_vec(1) + 1d0/4d0*EJ3xx * ZP*ZS
            Obs_scal(20)%Obs_vec(1)  =  Obs_scal(20)%Obs_vec(1) + 1d0/4d0*EJ3yy * ZP*ZS
            Obs_scal(21)%Obs_vec(1)  =  Obs_scal(21)%Obs_vec(1) + 1d0/4d0*EJ3zz * ZP*ZS
            

            do I = 1,Latt%N ! Sx,Sy,Sz
               do nc1 = 1, 2
                  select case(nc1)
                  case(1) ! Orbital A
                     I1_up = invlist(I,1)
                     I1_do = invlist(I,3)
                  case(2) ! Orbital B
                     I1_up = invlist(I,2)
                     I1_do = invlist(I,4)
                  end select
                  Esx=Esx+GRC(I1_up,I1_do,1)+GRC(I1_do,I1_up,1)
                  Esy=Esy+GRC(I1_up,I1_do,1)-GRC(I1_do,I1_up,1)
                  Esz=Esz+GRC(I1_up,I1_up,1)-GRC(I1_do,I1_do,1)
               end do
            end do   

            Obs_scal(9)%Obs_vec(1)  =  Obs_scal(9)%Obs_vec(1) + 1d0/2d0*(&
                 &  dble(Hx)*Esx&
                 & +dble(Hy)*cmplx(aimag(Esy),-real(Esy),kind(0.D0))&
                 & +dble(Hz)*Esz&
                 & ) * ZP*ZS
            
            Obs_scal(22)%Obs_vec(1)  =  Obs_scal(22)%Obs_vec(1) + 1d0/2d0*Esx * ZP*ZS
            Obs_scal(23)%Obs_vec(1)  =  Obs_scal(23)%Obs_vec(1) + 1d0/2d0*Esy * ZP*ZS
            Obs_scal(24)%Obs_vec(1)  =  Obs_scal(24)%Obs_vec(1) + 1d0/2d0*Esz * ZP*ZS
            
            
            do I=1,Latt%N
             I1=(I-1)*Norb
             do no_I=1,2

              do J=1,Latt%N
              J1=(J-1)*Norb
              do no_J=1,2
              imj = latt%imj(I,J)

              !xx 4*
              Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &+GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) )
              !yy (-4)*
              Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)-GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) )
              !zz 4*
              Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J+2,1)-GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I  , J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J  ,1)-GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J+2,1)+GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J+2,1) )

              !xy i4*
              Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) )
              !xz 4*
              Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(6)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J+2,1)-GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J+2,1)-GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J+2,1) )
              !yx i4*
              Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(7)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &+GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)-GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) )
              !yz i4*
              Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(8)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J+2,1)-GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J  ,1)-GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J+2,1)+GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I  , J1+no_J+2,1) )

              !zx 4*
              Obs_eq(9)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(9)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &+GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)-GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &-GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) )
              !zy i4*
              Obs_eq(10)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(10)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &-GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J  ,1)-GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I  , J1+no_J+2,1) &
              &-GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J+2,1)-GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J  ,1)+GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J+2,1) )

              !ni*nj
              Obs_eq(11)%Obs_Latt(imj,1,no_I,no_J)=Obs_eq(11)%Obs_Latt(imj,1,no_I,no_J)+Z*ZP*ZS*( &
              & GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I  , J1+no_J  ,1)*GR(I1+no_I  , J1+no_J  ,1) &
              &+GRC(I1+no_I  , I1+no_I  ,1)*GRC(J1+no_J+2, J1+no_J+2,1)+GRC(I1+no_I  , J1+no_J+2,1)*GR(I1+no_I  , J1+no_J+2,1) &
              &+GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J  , J1+no_J  ,1)+GRC(I1+no_I+2, J1+no_J  ,1)*GR(I1+no_I+2, J1+no_J  ,1) &
              &+GRC(I1+no_I+2, I1+no_I+2,1)*GRC(J1+no_J+2, J1+no_J+2,1)+GRC(I1+no_I+2, J1+no_J+2,1)*GR(I1+no_I+2, J1+no_J+2,1) )

              end do
              end do

              Obs_eq(11)%Obs_Latt0(no_I) =  Obs_eq(11)%Obs_Latt0(no_I) +Z*ZP*ZS*GRC(I1+no_I  , I1+no_I  ,1)
             end do
             end do

             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) =  Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) + &
                        &               Z * GRC(I1,J1,1) *  ZP*ZS
                ENDDO
             ENDDO

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
        Subroutine ObserT(NT,  GT0,G0T,G00,GTT, PHASE,  Mc_step_weight)

          Use Predefined_Obs

          Implicit none

          Integer         , INTENT(IN) :: NT
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL),G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
          Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
          Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight
          
          !Locals
          Complex (Kind=Kind(0.d0)) :: ZP, ZS
          ! Add local variables as needed
          Complex (Kind=Kind(0.d0)) :: Z
          Integer :: IMJ, I, J, I1, J1, no_I, no_J, Norb

          ZP = PHASE/Real(Phase, kind(0.D0))
          ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
          ZS = ZS * Mc_step_weight

          ! Compute observables
          If (NT == 0 ) then 
             DO I = 1,Size(Obs_tau,1)
                Obs_tau(I)%N = Obs_tau(I)%N + 1
                Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
             ENDDO
          endif

             Norb = Latt_unit%Norb
             Z =  cmplx(dble(N_SUN), 0.d0, kind(0.D0))
             do I=1,Latt%N
             I1=(I-1)*Norb
             do no_I=1,2

              do J=1,Latt%N
              J1=(J-1)*Norb
              do no_J=1,2
              imj = latt%imj(I,J)

              !xx 4*
              Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J  , J1+no_J+2,1)-G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J  , J1+no_J+2,1)-G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )
              !yy (-4)*
              Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &-GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J  , J1+no_J+2,1)+G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J+2, J1+no_J  ,1)+G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J  , J1+no_J+2,1)-G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )
              !zz 4*
              Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(4)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & (cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &           -G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &           +G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J+2,1) &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &           +G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &           -G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J+2,1) )


              !xy i4*
              Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(5)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &-GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J  , J1+no_J+2,1)+G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J  , J1+no_J+2,1)+G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )
              !xz 4* 
              Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(6)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              &-GTT(I1+no_I+2, I1+no_I  ,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &     -G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+GTT(I1+no_I+2, I1+no_I  ,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &     +G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &     -G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &     +G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )


              !yx i4*
              Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(7)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J+2, J1+no_J  ,1)-G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+GTT(I1+no_I+2, I1+no_I  ,1)*G00(J1+no_J  , J1+no_J+2,1)-G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J+2, J1+no_J  ,1)+G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*G00(J1+no_J  , J1+no_J+2,1)+G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )
              !yz i4*
              Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(8)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              &-GTT(I1+no_I+2, I1+no_I  ,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &     -G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+GTT(I1+no_I+2, I1+no_I  ,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &     +G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I+2, J1+no_J+2,1) &
              &+GTT(I1+no_I  , I1+no_I+2,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J  , J1+no_J  ,1))&
              &     +G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &-GTT(I1+no_I  , I1+no_I+2,1)*(cmplx(1.d0,0.d0,kind(0.d0))-G00(J1+no_J+2, J1+no_J+2,1))&
              &     -G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I  , J1+no_J+2,1) )

              !zx 4*
              Obs_tau(9)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(9)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*G00(J1+no_J+2, J1+no_J  ,1) &
              &      -G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*G00(J1+no_J  , J1+no_J+2,1) &
              &      -G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J+2,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*G00(J1+no_J+2, J1+no_J  ,1) &
              &      +G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*G00(J1+no_J  , J1+no_J+2,1) &
              &      +G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J+2,1) )
              !zy i4*
              Obs_tau(10)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(10)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*G00(J1+no_J+2, J1+no_J  ,1) &
              &      -G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*G00(J1+no_J  , J1+no_J+2,1) &
              &      +G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J+2,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*G00(J1+no_J+2, J1+no_J  ,1) &
              &      +G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &-(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*G00(J1+no_J  , J1+no_J+2,1) &
              &      -G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J+2,1) )

              !ni*nj
              Obs_tau(11)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(11)%Obs_Latt(imj,nt+1,no_I,no_J)+Z*ZP*ZS*( &
              & (cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*(1d0-G00(J1+no_J  , J1+no_J  ,1))&
              &           -G0T(J1+no_J  , I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J  ,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I  , I1+no_I  ,1))*(1d0-G00(J1+no_J+2, J1+no_J+2,1))&
              &           -G0T(J1+no_J+2, I1+no_I  ,1)*GT0(I1+no_I  , J1+no_J+2,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*(1d0-G00(J1+no_J  , J1+no_J  ,1))&
              &           -G0T(J1+no_J  , I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J  ,1) &
              &+(cmplx(1.d0,0.d0,kind(0.d0))-GTT(I1+no_I+2, I1+no_I+2,1))*(1d0-G00(J1+no_J+2, J1+no_J+2,1))&
              &           -G0T(J1+no_J+2, I1+no_I+2,1)*GT0(I1+no_I+2, J1+no_J+2,1) )

              end do
              end do


              Obs_tau(11)%Obs_Latt0(no_I) = Obs_tau(11)%Obs_Latt0(no_I) +Z*ZP*ZS* &
                     &         (cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1+no_I ,I1+no_I ,1))

             end do
             end do


             Do I1 = 1,Ndim
                I    = List(I1,1)
                no_I = List(I1,2)
                Do J1 = 1,Ndim
                   J    = List(J1,1)
                   no_J = List(J1,2)
                   imj = latt%imj(I,J)
                   ! Green
                   Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) =  Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J)  &
                        & +  Z * GT0(I1,J1,1) * ZP* ZS
                Enddo
             Enddo

        end Subroutine OBSERT
        
    end submodule ham_Kit_Heis_smod
