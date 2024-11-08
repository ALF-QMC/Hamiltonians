
   submodule (Hamiltonian_main) ham_AFM_N_S_mod

      Use Operator_mod
      Use Lattices_v3
      Use WaveFunction_mod
      ! Use Observables

      Implicit none

      integer, parameter :: dp=kind(0.d0)  ! double precision

      type, extends(ham_base) :: ham_AFM_N_S

      contains
         procedure, nopass :: Ham_Set                 => Ham_Set_AFM_N_S
         procedure, nopass :: Alloc_obs               => Alloc_obs_AFM_N_S
         procedure, nopass :: Obser                   => Obser_AFM_N_S
         procedure, nopass :: ObserT                  => ObserT_AFM_N_S
         !procedure, nopass :: Global_move             => Global_move_AFM_N_S
#ifdef HDF5
         procedure, nopass :: write_parameters_hdf5
#endif

      end type ham_AFM_N_S

      Integer                :: N_coord, Norb
      ! For orbital structure of Unit cell
      Integer, allocatable   :: List0(:,:), Invlist0(:,:), List2(:,:), Invlist2(:,:,:)
      Integer, allocatable   :: List_orb(:, :), Invlist_orb(:, :)

      Integer                :: N_Fam, L_intra, B_intra
      Integer, allocatable   :: checker_inter(:,:,:), checker_intra(:,:,:), L_FAM(:)

      !#PARAMETERS START# VAR_AFM_N_S
      Character (len=64) :: Lattice_type = 'Square'  ! Implemented: Square, N_leg_ladder and Open_square
      Integer :: L1 = 4                  ! Size of lattice in first dimension
      Integer :: L2 = 4                  ! Size of lattice in second dimension
      !Integer :: N_SUN = 2              ! SU(N) symmtry
      Integer :: N_Spin = 1              ! 2*S
      !logical :: Projector = .false.    ! Use projetive finite temperature method
      !logical :: Symm = .false.         ! Use symmetric Trotter decomposition
      real(dp) :: ham_J = 1.d0           ! Antiferromagnetic interaction strength
      real(dp) :: ham_U = 1.d0           ! On-site Hubbard term for freezing out charge fluctuations
      real(dp) :: ham_J2 = -1.d0         ! Projection parameter for going to fully symmetric representation
      real(dp) :: dtau = 0.1d0           ! Imaginary time step
      real(dp) :: beta = 5.d0            ! Reciprcal temperature
      real(dp) :: theta = 0.d0           ! Zero-temperature projection parameter
      real(dp) :: pinning_factor = 1.d0  ! Factor by which a single J-value will be multiplied to achieve phase pinning.
      integer :: pinning_x = 0           ! x coordinate of pinned bond
      integer :: pinning_y = 0           ! y coordinate of pinned bond
      !#PARAMETERS END#

#ifdef MPI
      Integer        :: Isize, Irank, irank_g, isize_g, igroup
      Integer        :: STATUS(MPI_STATUS_SIZE)
#endif

      Integer :: latt_id
      Type (Lattice),   target :: Latt
      Type (Unit_cell), target :: Latt_unit


   contains

   module subroutine Ham_Alloc_AFM_N_S()
     allocate(ham_AFM_N_S::ham)
   end subroutine Ham_Alloc_AFM_N_S

! Dynamically generated on compile time from parameters list.
! Supplies the subroutines read_parameters and write_parameters_hdf5.
#include "Hamiltonian_AFM_N_S_read_write_parameters.F90"

!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the Hamiltonian
!--------------------------------------------------------------------
      Subroutine Ham_Set_AFM_N_S()
         Implicit none
         Integer :: N_part
         
#ifdef MPI
         integer :: ierr
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERR)
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
         call MPI_Comm_rank(Group_Comm, irank_g, ierr)
         call MPI_Comm_size(Group_Comm, isize_g, ierr)
         igroup           = irank/isize_g
#endif


         ! From dynamically generated file "Hamiltonian_AFM_N_S_read_write_parameters.F90"
         call read_parameters()
         call write_info_AFM_N_S()

         Select case (Lattice_type)
         Case ("Square")
            latt_id = 1
         Case ("N_leg_ladder")
            latt_id = 2
         Case ("Open_square")
            latt_id = 3
         case default 
            Write(error_unit,*) "Ham_Set_AFM_N_S: Lattice '", Lattice_type, "' not yet implemented!"
            error stop
         end select

         !Set constant parameters
         N_FL = 1

         allocate(OP_T(1,1))
         call Op_make(Op_T(1,1),1)
         OP_T(1,1)%P(:) = 1
         OP_T(1,1)%O(:,:) = 0.d0
         OP_T(1,1)%g = 0.d0
         Call Op_set( Op_T(1,1) )

         call Ham_V()


         if (Projector) then
            N_part = Ndim/2
            Call Predefined_TrialWaveFunction(Lattice_type, L1, L2, N_part, N_FL, WF_L, WF_R, N_Spin)
            Thtrot = nint(theta/dtau)
#ifdef MPI
            If (Irank_g == 0 ) then
#endif
               print*, "gap: ", WF_L(1)%Degen, WF_R(1)%Degen
#ifdef MPI
            endif
#endif
         else
            Thtrot = 0
         endif
         Ltrot = nint(beta/dtau)
         Ltrot = Ltrot+2*Thtrot

         !call test_ham()

      end Subroutine Ham_Set_AFM_N_S

      Subroutine test_ham()
         Implicit none

         integer :: i, i_a, i_b
         integer :: unit_intra, unit_inter, unit
         real(dp) :: r_a(2), r_b(2)

         Open (newunit=unit_intra, file='test_intra.txt',status="unknown",position="append")
         Open (newunit=unit_inter, file='test_inter.txt',status="unknown",position="append")
         do i=1, size(op_v, 1)
            if (size(op_v(i, 1)%P) == 2) then
               i_a = op_v(i, 1)%P(1)
               i_b = op_v(i, 1)%P(2)
               if( (List2(i_a, 1) == List2(i_b, 1)) .and. (List2(i_a, 2) == List2(i_b, 2)) ) then
                  unit = unit_intra
               else
                  unit = unit_inter
               endif
               r_a = latt%list(List2(i_a, 1), 1) * latt%a1_p &
                   + latt%list(List2(i_a, 1), 2) * latt%a2_p
               r_a = r_a + latt_unit%Orb_pos_p(List2(i_a, 2), :)
               r_b = latt%list(List2(i_b, 1), 1) * latt%a1_p &
                   + latt%list(List2(i_b, 1), 2) * latt%a2_p
               r_b = r_b + latt_unit%Orb_pos_p(List2(i_b, 2), :)
               write(unit,*) op_v(i, 1)%P, op_v(i, 1)%g, r_a, r_b, List2(i_a, 3), List2(i_b, 3)
            endif
         enddo
         close(unit_intra)
         close(unit_inter)

      end Subroutine test_ham


      Subroutine write_info_AFM_N_S()
#ifdef MPI
         Use mpi
#endif
         Implicit none


         !Local variables
         integer            :: ierr
         Character (len=64) :: file_info

#ifdef MPI
         If (Irank_g == 0 ) then
#endif
#ifdef TEMPERING
            write(file_info,'(A,I0,A)') "Temp_",igroup,"/info"
#else
            file_info = "info"
#endif

            OPEN(Unit = 50,file=file_info,status="unknown",position="append")
            Write(50,*) '====================================='
            Write(50,*) 'Model is            : ', 'AFM_N_S'
            !Write(50,*) 'Lattice is          : ', Lattice_type
            !Write(50,*) '# of orbitals       : ', Ndim
            Write(50,*) 'Lattice_type        : ', Lattice_type
            Write(50,*) 'L1                  : ', L1
            Write(50,*) 'L2                  : ', L2
            Write(50,*) 'N_SUN               : ', N_SUN
            Write(50,*) 'N_Spin              : ', N_Spin
            Write(50,*) 'Symm                : ', Symm
            Write(50,*) 'Projector           : ', Projector
            Write(50,*) 'ham_J               : ', ham_J
            Write(50,*) 'ham_U               : ', ham_U
            Write(50,*) 'ham_J2              : ', ham_J2
            Write(50,*) 'dtau                : ', dtau
            Write(50,*) 'beta                : ', beta
            Write(50,*) 'theta               : ', theta
            Write(50,*) 'pinning_factor      : ', pinning_factor
            Write(50,*) 'pinning_x           : ', pinning_x
            Write(50,*) 'pinning_y           : ', pinning_y
            close(50)
#ifdef MPI
         endif
#endif
      end Subroutine write_info_AFM_N_S




!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Sets the interaction
!--------------------------------------------------------------------
      Subroutine Ham_V()

         Implicit none


         integer                                :: N_inter
         Complex (Kind=Kind(0.d0)), allocatable :: O_inter(:,:,:,:), g_inter(:,:), alpha_inter(:,:)
         integer, allocatable                   :: type_inter(:,:)

         integer                                :: N_intra
         Complex (Kind=Kind(0.d0)), allocatable :: O_intra(:,:,:,:), g_intra(:,:), alpha_intra(:,:)
         integer, allocatable                   :: type_intra(:,:)

         integer                                :: N_onsite
         Complex (Kind=Kind(0.d0)), allocatable :: O_onsite(:,:), g_onsite(:,:), alpha_onsite(:,:)
         integer, allocatable                   :: type_onsite(:,:)

         Integer, allocatable :: checker0(:,:,:), L_Fam0(:)
         
         Integer,                   allocatable :: I_overwrite(:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: overwrite_g_factor(:)

         N_inter = 2
         allocate( O_inter(2,2,N_inter,N_FL), g_inter(N_inter,N_FL) )
         allocate( alpha_inter(N_inter,N_FL), type_inter(N_inter,N_FL) )
         O_inter(:,:,:,:) = cmplx( 0.d0, 0.d0, kind(0.D0) )
         g_inter(1,:) = sqrt( cmplx( ham_J*dtau/8.d0, 0.d0, kind(0.D0) ) )
         g_inter(2,:) = sqrt( cmplx( ham_J*dtau/8.d0, 0.d0, kind(0.D0) ) )
         alpha_inter(:,:) = cmplx( 0.d0, 0.d0, kind(0.D0) )
         type_inter(:,:) = 2
         O_inter(1,2,1,1) = cmplx( 1.d0,  0.d0, kind(0.D0) )
         O_inter(2,1,1,1) = cmplx( 1.d0,  0.d0, kind(0.D0) )
         O_inter(1,2,2,1) = cmplx( 0.d0,  1.d0, kind(0.D0) )
         O_inter(2,1,2,1) = cmplx( 0.d0, -1.d0, kind(0.D0) )

         N_intra = 2
         allocate( O_intra(2,2,N_intra,N_FL), g_intra(N_intra,N_FL) )
         allocate( alpha_intra(N_intra,N_FL), type_intra(N_intra,N_FL) )
         O_intra(:,:,:,:) = cmplx( 0.d0, 0.d0, kind(0.D0) )
         g_intra(1,:) = sqrt( cmplx( ham_J2*dtau/8.d0, 0.d0, kind(0.D0) ) )
         g_intra(2,:) = sqrt( cmplx( ham_J2*dtau/8.d0, 0.d0, kind(0.D0) ) )
         alpha_intra(:,:) = cmplx( 0.d0, 0.d0, kind(0.D0) )
         type_intra(:,:) = 2
         O_intra(1,2,1,1) = cmplx( 1.d0,  0.d0, kind(0.D0) )
         O_intra(2,1,1,1) = cmplx( 1.d0,  0.d0, kind(0.D0) )
         O_intra(1,2,2,1) = cmplx( 0.d0,  1.d0, kind(0.D0) )
         O_intra(2,1,2,1) = cmplx( 0.d0, -1.d0, kind(0.D0) )

         N_onsite = 1
         allocate( O_onsite(N_onsite,N_FL), g_onsite(N_onsite,N_FL) )
         allocate( alpha_onsite(N_onsite,N_FL), type_onsite(N_onsite,N_FL) )
         O_onsite(1,:) = cmplx( 1.d0, 0.d0, kind(0.D0) )
         g_onsite(1,:) = sqrt(cmplx( -dtau*ham_U, 0.d0, kind(0.D0) ))
         alpha_onsite(1,:) = cmplx( -0.5d0, 0.d0, kind(0.D0) )
         type_onsite(1,:) = 2

         Call Predefined_Latt(Lattice_type, L1, L2, &
                              Ndim, List0, Invlist0, Latt, Latt_unit, &
                              List_orb, Invlist_orb)
         call Predefined_checkerboard(Lattice_type, Ndim, List0, Invlist0, Latt, Latt_unit, &
                                      L1, L2, List_orb, Invlist_orb, &
                                      checker0, L_Fam0)
         !call test_checkerboard(Lattice_type, Ndim, List0, Invlist0, Latt, Latt_unit, &
         !                       checker0, L_Fam0)
         call Checkeboard_Add_subfam(checker0, L_Fam0, latt%N, List0, Latt_unit%Norb, N_Spin, &
                                     List2, Invlist2, checker_inter, L_Fam, checker_intra)
         Ndim = Ndim*N_Spin
         N_Fam = size(checker_inter,3)
         B_intra = size(checker_intra,3)
         L_intra = size(checker_intra,2)

         call set_overwrite(I_overwrite, overwrite_g_factor)

         call NN_interaction( &
           Ndim, N_FL, checker_inter, L_Fam, checker_intra, Symm, &
           N_inter , O_inter , g_inter , alpha_inter , type_inter , &
           N_intra , O_intra , g_intra , alpha_intra , type_intra , &
           N_onsite, O_onsite, g_onsite, alpha_onsite, type_onsite, &
           I_overwrite, overwrite_g_factor, &
           OP_V )

         deallocate( O_inter, g_inter, alpha_inter, type_inter )
         deallocate( O_intra, g_intra, alpha_intra, type_intra )
         deallocate( O_onsite, g_onsite, alpha_onsite, type_onsite )
         deallocate( checker0, L_Fam0 )

      end Subroutine Ham_V


      subroutine test_checkerboard(Lattice_type, Ndim, List, Invlist, Latt, Latt_unit, &
                                 & checker, L_Fam)
         Implicit none

         Character (len=64), Intent(IN)           :: Lattice_type
         Integer, Intent(IN)                      :: Ndim
         Integer, Intent(IN), Dimension(:,:)      :: List, Invlist
         Type(Lattice), Intent(in)                :: Latt
         Type(Unit_cell), Intent(in)              :: Latt_unit
         
         Integer, Intent(in) :: checker(:,:,:)
         integer, Intent(in) :: L_Fam(:)

         integer :: a, b, I1, J1, no1, I2, J2, no2, unit
         real(dp) :: r1(2), r2(2)

         Open (newunit=unit, file='test_checker.txt', status="unknown", position="append")
         do a=1, size(checker,3)
            do b=1, L_FAM(a)
               I1 = checker(1, b, a)
               J1 = List(I1, 1)
               no1 = List(I1, 2)
               r1 = latt%list(J1, 1) * latt%a1_p + latt%list(J1, 2) * latt%a2_p + Latt_unit%Orb_pos_p(no1, :)

               I2 = checker(2, b, a)
               J2 = List(I2, 1)
               no2 = List(I2, 2)
               r2 = latt%list(J2, 1) * latt%a1_p + latt%list(J2, 2) * latt%a2_p + Latt_unit%Orb_pos_p(no2, :)
               write(unit, *) r1, r2
            enddo
         enddo
         close(unit)
      end subroutine test_checkerboard


      subroutine set_overwrite(I_overwrite, overwrite_g_factor)
         implicit none
         Integer,                   allocatable, intent(inout) :: I_overwrite(:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable, intent(inout) :: overwrite_g_factor(:)

         Integer :: I1, I2, no1, no2, nc
         Integer :: J1, J2, a, b, x1, x2, y1, y2, unit
         Integer :: norb_pin1, norb_pin2


         
         allocate(I_overwrite(2, N_Spin**2, N_FL), overwrite_g_factor(N_FL))
         if (Lattice_type == "Open_square") then
            !norb_pin1 = Invlist_orb(L1/2, 1)
            !norb_pin2 = Invlist_orb(L1/2+1, 1)
            norb_pin1 = Invlist_orb(pinning_x+1, pinning_y+1)
            norb_pin2 = Invlist_orb(pinning_x+2, pinning_y+1)
            print*, "norb_pins", norb_pin1, norb_pin2
            print*, "coord1", latt_unit%Orb_pos_p(norb_pin1, 1), latt_unit%Orb_pos_p(norb_pin1, 2)
            print*, "coord2", latt_unit%Orb_pos_p(norb_pin2, 1), latt_unit%Orb_pos_p(norb_pin2, 2)
            nc = 0
            do no1 = 1, N_Spin
               do no2 = 1, N_Spin
                  nc = nc+1
                  I_overwrite(1, nc, :) = Invlist2(1, norb_pin1, no1)
                  I_overwrite(2, nc, :) = Invlist2(1, norb_pin2, no2)
               enddo
            enddo
         else if (Lattice_type == "N_leg_ladder") then
            I1 = latt%Invlist(pinning_x, 0)
            I2 = latt%Invlist(pinning_x+1, 0)
            nc = 0
            do no1 = 1, N_Spin
               do no2 = 1, N_Spin
                  nc = nc+1
                  I_overwrite(1, nc, :) = Invlist2(I1, pinning_y+1, no1)
                  I_overwrite(2, nc, :) = Invlist2(I2, pinning_y+1, no2)
               enddo
            enddo
         else if (Lattice_type == "Square") then
            I1 = latt%Invlist(pinning_x, pinning_y)
            I2 = latt%Invlist(pinning_x+1, pinning_y)
            nc = 0
            do no1 = 1, N_Spin
               do no2 = 1, N_Spin
                  nc = nc+1
                  I_overwrite(1, nc, :) = Invlist2(I1, 1, no1)
                  I_overwrite(2, nc, :) = Invlist2(I2, 1, no2)
               enddo
            enddo
         else
            write(error_unit, *) "Pinning not defined for lattice " // Lattice_type
            error stop
         endif
         overwrite_g_factor(:) = sqrt(pinning_factor)
         
         ! allocate(I_overwrite(2, N_Spin**2*L2, N_FL), overwrite_g_factor(N_FL))
         ! Open (newunit=unit, file='test_bonds1.txt',status="unknown",position="append")
         ! nc = 0
         ! do a=1, N_Fam
         !    do b=1, L_FAM(a)
         !       I1 = checker_inter(1, b, a)
         !       I2 = checker_inter(2, b, a)
         !       J1 = List2(I1, 1)
         !       x1 = latt%list(J1, 1)
         !       y1 = latt%list(J1, 2)
         !       J2 = List2(I2, 1)
         !       x2 = latt%list(J2, 1)
         !       y2 = latt%list(J2, 2)
         !       if( x1 == x2 .and. ((y1 == 0 .and. y2 == 1) .or. (y1 == 1 .and. y2 == 0)) ) then
         !          nc = nc+1
         !          I_overwrite(1, nc, :) = I1
         !          I_overwrite(2, nc, :) = I2
         !          write(unit,*) I1, I2
         !       endif
         !    enddo
         ! enddo
         ! close(unit)
         ! overwrite_g_factor(:) = 1E-6

      end subroutine set_overwrite
      
      
      Subroutine NN_interaction( &
        Ndim, N_FL, checker_inter, L_Fam, checker_intra, Symm, &
        N_inter , O_inter , g_inter , alpha_inter , type_inter , &
        N_intra , O_intra , g_intra , alpha_intra , type_intra , &
        N_onsite, O_onsite, g_onsite, alpha_onsite, type_onsite, &
        I_overwrite, overwrite_g_factor, &
        OP_V )
        
        Implicit none
        
        Integer,                   Intent(in) :: Ndim, N_FL, checker_inter(:,:,:), L_Fam(:), checker_intra(:,:,:)
        Logical,                   Intent(in) :: Symm
        Integer,                   Intent(in) :: N_inter, N_intra, N_onsite
        Complex (Kind=Kind(0.d0)), Intent(in) :: O_inter(2,2,N_inter,N_FL), g_inter(N_inter,N_FL), alpha_inter(N_inter,N_FL)
        Complex (Kind=Kind(0.d0)), Intent(in) :: O_intra(2,2,N_inter,N_FL), g_intra(N_inter,N_FL), alpha_intra(N_inter,N_FL)
        Complex (Kind=Kind(0.d0)), Intent(in) :: O_onsite(N_onsite,N_FL), g_onsite(N_onsite,N_FL), alpha_onsite(N_onsite,N_FL)
        Integer,                   Intent(in) :: type_inter(N_inter,N_FL), type_intra(N_inter,N_FL), type_onsite(N_onsite,N_FL)
        
        Integer,                   Intent(in) :: I_overwrite(:, :, :)
        Complex (Kind=Kind(0.d0)), Intent(in) :: overwrite_g_factor(N_FL)
        
        Type(Operator), Intent(Out), allocatable :: Op_V(:,:) 

        !Local
        Integer :: I, nf, nc, nc1, nc2, I_term
        Integer :: N_Fam0,  N_int
        Integer :: B_intra, L_intra, N_int2
        Complex (Kind=Kind(0.d0)) :: f

        !integer :: n_test
        !n_test = 0
        
        ! Integer :: I_Op
        !print*, "I_overwrite", I_overwrite
        
        
        N_Fam0 = size(checker_inter,3)
        B_intra= size(checker_intra,3)
        L_intra= size(checker_intra,2)
                              
        If (Symm) then
          N_int = 0
          do nc = 1, N_Fam0
            N_int = N_int +  2 * N_inter * L_FAM(nc)
          enddo
          N_int2 = 2 * N_intra * B_intra * L_intra
        else
          N_int = 0
          do nc = 1, N_Fam0
            N_int = N_int +  N_inter * L_FAM(nc)
          enddo
          N_int2 = N_intra * B_intra * L_intra
        endif
        
        Allocate(Op_V(N_int + N_int2 + Ndim*N_onsite,N_FL))
        
        do nf = 1,N_FL
          nc = 0
          
          ! Set onsite interactions
          do I = 1,Ndim
            do I_term = 1, N_onsite
              nc = nc + 1
              call Op_make(Op_V(nc,nf),1)
              Op_V(nc,nf)%P(1)   = I
              Op_V(nc,nf)%O(1,1) = O_onsite(I_term,nf)
              Op_V(nc,nf)%g      = g_onsite(I_term,nf)
              Op_V(nc,nf)%alpha  = alpha_onsite(I_term,nf)
              Op_V(nc,nf)%type   = type_onsite(I_term,nf)
              Call Op_set( Op_V(nc,nf) )
            enddo
          enddo
          
          ! Set intrasite interactions
          if (Symm) then
            f = 1.0d0
          else
            f = sqrt(2.0d0)
          endif
          do nc1 = 1,B_intra
            do I_term = 1, N_intra
              do nc2 = 1,L_intra
                nc = nc + 1
                call Op_make(Op_V(nc,nf),2)
                Op_V(nc,nf)%P = checker_intra(:, nc2, nc1)
                Op_V(nc,nf)%O = O_intra(:,:,I_term,nf)
                Op_V(nc,nf)%g      = f*g_intra(I_term,nf)
                Op_V(nc,nf)%alpha  = alpha_intra(I_term,nf)
                Op_V(nc,nf)%type   = type_intra(I_term,nf)
                Call Op_set( Op_V(nc,nf) )
              enddo
            enddo
          enddo
          
          ! Set intersite interactions
          if (Symm) then
            f = sqrt(0.5d0)
          else
            f = 1.0d0
          endif
          do nc1 = 1,N_Fam0
            do I_term = 1, N_inter
              do nc2 = 1, L_FAM(nc1)
                nc = nc + 1
                call Op_make(Op_V(nc,nf),2)
                Op_V(nc,nf)%P = checker_inter(:, nc2, nc1)
                Op_V(nc,nf)%O = O_inter(:,:,I_term,nf)
                Op_V(nc,nf)%g      = f*g_inter(I_term,nf)
                do I=1, size(I_overwrite, 2)
                  if ((I_overwrite(1, I, nf) == checker_inter(1, nc2, nc1) .and. &
                       I_overwrite(2, I, nf) == checker_inter(2, nc2, nc1)) .or. &
                      (I_overwrite(1, I, nf) == checker_inter(2, nc2, nc1) .and. &
                       I_overwrite(2, I, nf) == checker_inter(1, nc2, nc1))) then
                    Op_V(nc,nf)%g = Op_V(nc,nf)%g * overwrite_g_factor(nf)
                    !n_test = n_test+1
                    !print*, "test", n_test, nc1, nc2
                  endif
                enddo
                Op_V(nc,nf)%alpha  = alpha_inter(I_term,nf)
                Op_V(nc,nf)%type   = type_inter(I_term,nf)
                Call Op_set( Op_V(nc,nf) )
               !  Op_V(nc,nf)%I_Fam = nc1
              enddo
            enddo
          enddo
          
          ! Mirror interactions in case of symmetrization
          if (Symm) then
            nc2 = nc
            do nc1 = 1, N_int/2 + N_int2/2
              nc = nc + 1
              call Op_make(Op_V(nc,nf), 2)
              Op_V(nc, nf)%P = Op_V(nc2, nf)%P
              Op_V(nc, nf)%O = Op_V(nc2, nf)%O
              Op_V(nc,nf)%g = Op_V(nc2, nf)%g
              Op_V(nc,nf)%alpha = Op_V(nc2, nf)%alpha
              Op_V(nc,nf)%type = Op_V(nc2, nf)%type
              Call Op_set( Op_V(nc,nf) )
              nc2 = nc2 - 1
            enddo
          endif
          
          if ( nc /= size(Op_V,1) ) then
            Write(*,*) "Error in NN_interaction", nc, size(Op_V,1)
            stop
          endif
        enddo
        
      End Subroutine NN_interaction



      Subroutine Obser_Latt_make_norb(Obs, Nt, Norb, Filename, Latt, Channel, dtau)
         Implicit none
         class(Obser_Latt), Intent(INOUT)      :: Obs
         Integer,           Intent(IN)         :: Nt
         Integer,           Intent(IN)         :: Norb
         Character(len=64), Intent(IN)         :: Filename
         Type(Lattice),     Intent(IN), target :: Latt
         Character(len=2),  Intent(IN)         :: Channel
         Real(Kind=Kind(0.d0)),  Intent(IN)    :: dtau
        
         allocate(Obs%Latt_unit)
         Obs%Latt_unit%N_coord = 0
         Obs%Latt_unit%Norb = Norb
         allocate(Obs%Latt_unit%Orb_pos_p(Norb, 2))
         Obs%Latt_unit%Orb_pos_p = 0.d0
         
         Allocate (Obs%Obs_Latt(Latt%N, Nt, Obs%Latt_unit%Norb, Obs%Latt_unit%Norb))
         Allocate (Obs%Obs_Latt0(Obs%Latt_unit%Norb))
         Obs%File_Latt = Filename
         Obs%Latt => Latt
         Obs%Channel = Channel
         Obs%dtau = dtau
      end subroutine Obser_Latt_make_norb


!--------------------------------------------------------------------
!> @author
!> ALF Collaboration
!>
!> @brief
!> Specifiy the equal time and time displaced observables
!> @details
!--------------------------------------------------------------------
      Subroutine  Alloc_obs_AFM_N_S(Ltau)

         Implicit none
         !>  Ltau=1 if time displaced correlations are considered.
         Integer, Intent(In) :: Ltau
         Integer    ::  i, N, Nt, No
         Character (len=64) :: Filename
         Character (len=2)  :: Channel


         ! Scalar observables
         Allocate ( Obs_scal(4) )
         Do I = 1,Size(Obs_scal,1)
            select case (I)
            case (1)
               N = 1;   Filename ="n_test"
            case (2)
               N = 2;   Filename ="Ener"
            case (3)
               N = 3;   Filename ="S_squared"
            case (4)
               N = 2*Latt_unit%Norb;   Filename ="bonds_back"
            case default
               Write(6,*) ' Error in Alloc_obs '
            end select
            Call Obser_vec_make(Obs_scal(I),N,Filename)
         enddo


         ! Equal time correlators
         Allocate ( Obs_eq(11) )
         Do I = 1,Size(Obs_eq,1)
            select case (I)
            case (1)
               No = Latt_unit%Norb;  Filename ="Green"
            case (2)
               No = Latt_unit%Norb;  Filename ="Den"
            case (3)
               No = Latt_unit%Norb;  Filename ="Spin"
            case (4)
               No = Latt_unit%Norb;  Filename ="Spin2"
            case (5)
               No = Latt_unit%Norb;  Filename ="Spin3"
            case (6)
               No = 2*Latt_unit%Norb;     Filename ="Dimer"
            case (7)
               No = 2*Latt_unit%Norb;     Filename ="Dimer2"
            case (8)
               No = 2*Latt_unit%Norb;     Filename ="Dimer3"
            case (9)
               No = 2*Latt_unit%Norb;     Filename ="Dimer_Z"
            case (10)
               No = 2*Latt_unit%Norb;     Filename ="bonds"
            case (11)
               No = 2*Latt_unit%Norb;     Filename ="bonds_"
            case default
               Write(6,*) ' Error in Alloc_obs '
            end select
            Nt = 1
            Channel = '--'
            call Obser_Latt_make_norb(Obs_eq(I), Nt, No, Filename, Latt, Channel, dtau)
         enddo

         If (Ltau == 1) then
            ! Time-displaced correlators
            Allocate ( Obs_tau(3) )
            Do I = 1,Size(Obs_tau,1)
               select case (I)
               case (1)
                  Channel = 'P';  Filename ="Green"
               case (2)
                  Channel = 'PH'; Filename ="Den"
               case (3)
                  Channel = 'PH'; Filename ="Spin"
               case default
                  Write(6,*) ' Error in Alloc_obs '
               end select
               Nt = Ltrot+1-2*Thtrot
               call Obser_Latt_make(Obs_tau(I), Nt, Filename, Latt, Latt_unit, Channel, dtau)
            enddo
         endif
      End Subroutine Alloc_obs_AFM_N_S



      pure function spin(GR, GRC, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: spin

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         spin = NSUN*cmplx(0.5d0, 0.d0, kind(0.D0)) * &
         & ( -GR(j,j,nf)*GRC(i,i,nf) -NSUN*GRC(i,j,nf)*GRC(j,i,nf) + cmplx(0.25d0, 0.d0, kind(0.D0)) )
      end function spin

      pure function spin_tau(GT0, G0T, G00, GTT, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: spin_tau

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(:,:,:),G0T(:,:,:),G00(:,:,:),GTT(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         spin_tau = NSUN*cmplx(0.5d0, 0.d0, kind(0.D0)) * &
         & ( -G00(j,j,nf)*( cmplx(1.d0, 0.d0, kind(0.D0))-GTT(i,i,nf) ) &
             -NSUN*GT0(i,j,nf)*G0T(j,i,nf) + cmplx(0.25d0, 0.d0, kind(0.D0)) )
      end function spin_tau

      pure function spin2(GR, GRC, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: spin2

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         spin2 = ((-1 + NSUN**2)*GR(i,j,nf)*GRC(i,j,nf))/2.
      end function spin2

      pure function dimer2(GR, GRC, NSUN, i, j, k, l)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: dimer2

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j, k, l

         Integer, parameter :: nf = 1

         dimer2 = ((-1 + NSUN**2)*(GR(i,k,nf)*(GR(j,l,nf)*GRC(i,l,nf)*GRC(j,k,nf) + &
          &         (-((-1 + NSUN**2)*GR(k,l,nf)*GRC(i,j,nf)) + NSUN*GR(j,l,nf)*GRC(i,k,nf))*GRC(j,l,nf)) + &
          &      (-1 + NSUN**2)*GR(i,j,nf)* &
          &       (GR(j,k,nf)*GR(k,l,nf)*GRC(i,l,nf) + &
          &         (NSUN*GR(k,l,nf)*GRC(i,j,nf) - GR(j,l,nf)*GRC(i,k,nf))*GRC(k,l,nf)) + &
          &       GR(i,l,nf)*(GR(j,k,nf)*(NSUN*GRC(i,l,nf)*GRC(j,k,nf) + &
          &            GRC(i,k,nf)*GRC(j,l,nf)) + &
          &         (-1 + NSUN**2)*GRC(i,j,nf)*GRC(j,k,nf)*GRC(k,l,nf))))/(4.*NSUN)
      end function dimer2

      pure function spin_z(GR, GRC, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: spin_z

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:)
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         spin_z = 2*GR(i,j,nf)*GRC(i,j,nf)
      end function spin_z

      pure function dimer_z(GR, GRC, i, j, k, l)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: dimer_z

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:)
         Integer,                   INTENT(IN) :: i, j, k, l

         Integer, parameter :: nf = 1

         dimer_z = 2*(-(GR(i,k,nf)*(GR(j,l,nf)*GRC(i,l,nf)*GRC(j,k,nf) + &
          &         (GR(k,l,nf)*GRC(i,j,nf) - 2*GR(j,l,nf)*GRC(i,k,nf))*GRC(j,l,nf))) + &
          &    GR(i,j,nf)*(GR(j,k,nf)*GR(k,l,nf)*GRC(i,l,nf) + &
          &       (2*GR(k,l,nf)*GRC(i,j,nf) - GR(j,l,nf)*GRC(i,k,nf))*GRC(k,l,nf)) + &
          &    GR(i,l,nf)*(GR(j,k,nf)*(2*GRC(i,l,nf)*GRC(j,k,nf) - &
          &          GRC(i,k,nf)*GRC(j,l,nf)) + GRC(i,j,nf)*GRC(j,k,nf)*GRC(k,l,nf)))
      end function dimer_z

      pure function spin3(GR, GRC, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: spin3

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         spin3 = cmplx(0.5d0,0.d0,kind(0.D0))*NSUN* &
            ( NSUN*GR(i,j,nf)*GRC(i,j,nf) + GRC(i,i,nf)*GRC(j,j,nf) -cmplx(0.25d0,0.d0,kind(0.D0)) )
      end function spin3

      pure function dimer3(GR, GRC, NSUN, i, j, k, l)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: dimer3

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j, k, l

         Integer, parameter :: nf = 1

         dimer3 = cmplx(0.25d0,0.d0,kind(0.D0))*NSUN* &
          &   (NSUN**2*GR(i,j,nf)*GR(j,k,nf)*GR(k,l,nf)*GRC(i,l,nf) + &
          &    NSUN*GR(i,k,nf)*GR(k,l,nf)*GRC(i,l,nf)*GRC(j,j,nf) + &
          &    NSUN*GR(i,l,nf)*GR(j,k,nf)*GRC(i,l,nf)*GRC(j,k,nf) - &
          &    GR(i,k,nf)*GR(j,l,nf)*GRC(i,l,nf)*GRC(j,k,nf) + &
          &    NSUN*GR(j,k,nf)*GR(k,l,nf)*GRC(i,i,nf)*GRC(j,l,nf) - &
          &    NSUN**2*GR(i,k,nf)*GR(k,l,nf)*GRC(i,j,nf)*GRC(j,l,nf) - &
          &    GR(i,l,nf)*GR(j,k,nf)*GRC(i,k,nf)*GRC(j,l,nf) + &
          &    NSUN*GR(i,k,nf)*GR(j,l,nf)*GRC(i,k,nf)*GRC(j,l,nf) + &
          &    NSUN*GR(i,j,nf)*GR(j,l,nf)*GRC(i,l,nf)*GRC(k,k,nf) + &
          &    GR(i,l,nf)*GRC(i,l,nf)*GRC(j,j,nf)*GRC(k,k,nf) + &
          &    GR(j,l,nf)*GRC(i,i,nf)*GRC(j,l,nf)*GRC(k,k,nf) - &
          &    NSUN*GR(i,l,nf)*GRC(i,j,nf)*GRC(j,l,nf)*GRC(k,k,nf) + &
          &    NSUN**3*GR(i,j,nf)*GR(k,l,nf)*GRC(i,j,nf)*GRC(k,l,nf) - &
          &    NSUN**2*GR(i,j,nf)*GR(j,l,nf)*GRC(i,k,nf)*GRC(k,l,nf) + &
          &    NSUN**2*GR(k,l,nf)*GRC(i,i,nf)*GRC(j,j,nf)*GRC(k,l,nf) - &
          &    NSUN*GR(i,l,nf)*GRC(i,k,nf)*GRC(j,j,nf)*GRC(k,l,nf) - &
          &    NSUN*GR(j,l,nf)*GRC(i,i,nf)*GRC(j,k,nf)*GRC(k,l,nf) + &
          &    NSUN**2*GR(i,l,nf)*GRC(i,j,nf)*GRC(j,k,nf)*GRC(k,l,nf) + &
          &    NSUN*GR(i,j,nf)*GR(j,k,nf)*GRC(i,k,nf)*GRC(l,l,nf) + &
          &    GR(i,k,nf)*GRC(i,k,nf)*GRC(j,j,nf)*GRC(l,l,nf) + &
          &    GR(j,k,nf)*GRC(i,i,nf)*GRC(j,k,nf)*GRC(l,l,nf) - &
          &    NSUN*GR(i,k,nf)*GRC(i,j,nf)*GRC(j,k,nf)*GRC(l,l,nf) + &
          &    NSUN**2*GR(i,j,nf)*GRC(i,j,nf)*GRC(k,k,nf)*GRC(l,l,nf) + &
          &    NSUN*GRC(i,i,nf)*GRC(j,j,nf)*GRC(k,k,nf)*GRC(l,l,nf))
      end function dimer3

      pure function dimer3_0(GR, GRC, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: dimer3_0

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         dimer3_0 = cmplx(0.5d0,0.d0,kind(0.D0))*NSUN* (NSUN*GR(i,j,nf)*GRC(i,j,nf) + GRC(i,i,nf)*GRC(j,j,nf))
      end function dimer3_0

      pure function DDDD(GR, GRC, NSUN, i, j, k, l)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: DDDD

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j, k, l

         Integer, parameter :: nf = 1

         DDDD = NSUN*(NSUN*GR(j,j,nf)*GR(l,l,nf)*GRC(i,i,nf)*GRC(k,k,nf) - &
          &    GR(i,l,nf)*GR(j,j,nf)*GRC(i,l,nf)*GRC(k,k,nf) + &
          &    NSUN**2*GR(l,l,nf)*GRC(i,j,nf)*GRC(j,i,nf)*GRC(k,k,nf) - &
          &    NSUN*GR(j,l,nf)*GRC(i,l,nf)*GRC(j,i,nf)*GRC(k,k,nf) + &
          &    GR(j,l,nf)*GRC(i,i,nf)*GRC(j,l,nf)*GRC(k,k,nf) - &
          &    NSUN*GR(i,l,nf)*GRC(i,j,nf)*GRC(j,l,nf)*GRC(k,k,nf) + &
          &    NSUN*GR(i,l,nf)*GR(j,j,nf)*GRC(i,k,nf)*GRC(k,l,nf) + &
          &    NSUN**2*GR(j,l,nf)*GRC(i,k,nf)*GRC(j,i,nf)*GRC(k,l,nf) - &
          &    NSUN*GR(j,l,nf)*GRC(i,i,nf)*GRC(j,k,nf)*GRC(k,l,nf) + &
          &    NSUN**2*GR(i,l,nf)*GRC(i,j,nf)*GRC(j,k,nf)*GRC(k,l,nf) + &
          &    NSUN**2*(GR(j,j,nf)*GRC(i,i,nf) + NSUN*GRC(i,j,nf)*GRC(j,i,nf))*&
          &     GRC(k,l,nf)*GRC(l,k,nf) + &
          &    GR(j,k,nf)*(GR(l,l,nf)*(NSUN*GRC(i,k,nf)*GRC(j,i,nf) - &
          &          GRC(i,i,nf)*GRC(j,k,nf)) + &
          &       GR(i,l,nf)*(NSUN*GRC(i,l,nf)*GRC(j,k,nf) - &
          &          GRC(i,k,nf)*GRC(j,l,nf)) + &
          &       NSUN*(NSUN*GRC(i,l,nf)*GRC(j,i,nf) - GRC(i,i,nf)*GRC(j,l,nf))*&
          &        GRC(l,k,nf)) + GR(i,k,nf)*&
          &     (GR(j,j,nf)*GR(l,l,nf)*GRC(i,k,nf) + &
          &       NSUN*GR(l,l,nf)*GRC(i,j,nf)*GRC(j,k,nf) - &
          &       GR(j,l,nf)*GRC(i,l,nf)*GRC(j,k,nf) + &
          &       NSUN*GR(j,l,nf)*GRC(i,k,nf)*GRC(j,l,nf) + &
          &       NSUN*(GR(j,j,nf)*GRC(i,l,nf) + NSUN*GRC(i,j,nf)*GRC(j,l,nf))*&
          &        GRC(l,k,nf)))
      end function DDDD

      pure function DD(GR, GRC, NSUN, i, j)
         Implicit none
         Complex (Kind=Kind(0.d0)) :: DD

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(:,:,:), GRC(:,:,:), NSUN
         Integer,                   INTENT(IN) :: i, j

         Integer, parameter :: nf = 1

         DD = NSUN*(GR(j,j,nf)*GRC(i,i,nf) + NSUN*GRC(i,j,nf)*GRC(j,i,nf))
      end function DD


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
      subroutine Obser_AFM_N_S(GR,Phase,Ntau, Mc_step_weight)

         Implicit none

         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GR(Ndim,Ndim,N_FL)
         Complex (Kind=Kind(0.d0)), Intent(IN) :: PHASE
         Integer, INTENT(IN)          :: Ntau
         Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

         !Local
         Complex (Kind=Kind(0.d0)) :: GRC(Ndim,Ndim,N_FL), ZP, ZS
         Integer :: I, J, nf, I1, J1, imj, no_I, no_J, no2_I, no2_J, nc1, nc2
         Complex (Kind=Kind(0.d0)) :: Z, Z2, Z3, ZN, ZS0_sq
         Complex (Kind=Kind(0.d0)), allocatable :: Z_back(:)

         Integer :: I2, J2, n1, n2, no2_I2, no2_J2, nI, nJ
         Integer :: I1x, I1y, I2x, I2y, J1x, J1y, J2x, J2y
         Integer :: no_I2, no_J2

         ZP = PHASE/Real(Phase, kind(0.D0))
         ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))

         ZN = cmplx(dble(N_SUN), 0.d0, kind(0.D0))

         ZS0_sq = cmplx( ( dble(N_SUN)**2 + dble(N_SUN) )/8.d0 , 0.d0, kind(0.D0))


         Do nf = 1,N_FL
            Do I = 1,Ndim
               Do J = 1,Ndim
                  GRC(I, J, nf) = -GR(J, I, nf)
               Enddo
               GRC(I, I, nf) = 1.D0 + GRC(I, I, nf)
            Enddo
         Enddo
         ! GRC(i,j,nf) = < c^{dagger}_{j,nf } c_{j,nf } >

         ! Compute scalar observables.
         !Do I = 1,Size(Obs_scal,1)
         !   Obs_scal(I)%N         =  Obs_scal(I)%N + 1
         !   Obs_scal(I)%Ave_sign  =  Obs_scal(I)%Ave_sign + Real(ZS,kind(0.d0))
         !Enddo

         !n_test
         Z = cmplx(0.d0, 0.d0, kind(0.D0))
         Do I = 1,Ndim
            Z = Z + ZN * ( GRC(I, I, 1) - cmplx(0.5d0, 0.d0, kind(0.D0)) )**2 + GRC(I, I, 1)*GR(I,I,1)
         Enddo
         !Obs_scal(1)%Obs_vec(1) = Obs_scal(1)%Obs_vec(1) + Z*ZN*ZP*ZS/cmplx(dble(Ndim), 0.d0, kind(0.D0))
         Z = Z*ZN/cmplx(dble(Ndim), 0.d0, kind(0.D0))
         call Obs_scal(1)%measure( (/Z/), Phase )

         !Energy
         Z = cmplx(0.d0, 0.d0, kind(0.D0))
         Z2 = cmplx(0.d0, 0.d0, kind(0.D0))
         do nc1 = 1,N_Fam
            do nc2 = 1,L_FAM(nc1)
               I = checker_inter(1, nc2, nc1)
               J = checker_inter(2, nc2, nc1)
               Z = Z + spin(GR, GRC, ZN, I, J)
               Z2 = Z2 + spin2(GR, GRC, ZN, I, J)
            enddo
         enddo
         !Obs_scal(2)%Obs_vec(1) = Obs_scal(2)%Obs_vec(1) + Z*ZP*ZS*cmplx(ham_J, 0.d0, kind(0.D0))
         !Obs_scal(2)%Obs_vec(2) = Obs_scal(2)%Obs_vec(2) + Z2*ZP*ZS*cmplx(ham_J, 0.d0, kind(0.D0))
         Z  = Z *cmplx(ham_J, 0.d0, kind(0.D0))
         Z2 = Z2*cmplx(ham_J, 0.d0, kind(0.D0))
         call Obs_scal(2)%measure( (/Z, Z2/), Phase )

         !S squared
         Z = cmplx(0.d0, 0.d0, kind(0.D0))
         Z2 = cmplx(0.d0, 0.d0, kind(0.D0))
         Z3 = cmplx(0.d0, 0.d0, kind(0.D0))
         do nc1 = 1,B_intra
            do nc2 = 1,L_intra
               I = checker_intra(1, nc2, nc1)
               J = checker_intra(2, nc2, nc1)
               Z = Z + spin(GR, GRC, ZN, I, J)
               Z2 = Z2 + spin2(GR, GRC, ZN, I, J)
               Z3 = Z3 + spin3(GR, GRC, ZN, I, J)
            enddo
         enddo

         Z  = Z *cmplx(2.d0/dble(Latt%N*Latt_unit%Norb), 0.d0, kind(0.D0)) &
             + ZS0_sq*cmplx(dble(N_Spin), 0.d0, kind(0.D0)) !Factor 2 included
         Z2 = Z2*cmplx(2.d0/dble(Latt%N*Latt_unit%Norb), 0.d0, kind(0.D0)) &
             + ZS0_sq*cmplx(dble(N_Spin), 0.d0, kind(0.D0)) !Factor 2 included
         Z3 = Z3*cmplx(2.d0/dble(Latt%N*Latt_unit%Norb), 0.d0, kind(0.D0)) &
             + ZS0_sq*cmplx(dble(N_Spin), 0.d0, kind(0.D0)) !Factor 2 included
         call Obs_scal(3)%measure( (/Z, Z2, Z3/), Phase )


         ! Compute correlation functions
         DO I = 1,Size(Obs_eq,1)
            Obs_eq(I)%N        = Obs_eq(I)%N + 1
            Obs_eq(I)%Ave_sign = Obs_eq(I)%Ave_sign + real(ZS,kind(0.d0))
         ENDDO
         allocate(Z_back(2*Latt_unit%Norb))
         Z_back = cmplx(0.d0, 0.d0, kind(0.D0))
         Do I1 = 1,Ndim
            I     = List2(I1,1)
            no_I  = List2(I1,2)
            no2_I = List2(I1,3)
            Do J1 = 1,Ndim
               J     = List2(J1,1)
               no_J  = List2(J1,2)
               no2_J = List2(J1,3)
               imj = latt%imj(I,J)
               ! Green
               Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(1)%Obs_Latt(imj,1,no_I,no_J) &
                  & + ZN * GRC(I1,J1,1) * ZP*ZS

               ! Den
               Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(2)%Obs_Latt(imj,1,no_I,no_J) &
                  & + ( ZN*GRC(I1,I1,1) * GRC(J1,J1,1) &
                  &     +  GRC(I1,J1,1) * GR(I1,J1,1 )   ) * ZN * ZP*ZS

               ! Spin
               if ( I1 == J1 ) then
                  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) &
                  & + ZP*ZS * ZS0_sq
               else
                  Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(3)%Obs_Latt(imj,1,no_I,no_J) &
                  & + ZP*ZS * spin(GR, GRC, ZN, I1, J1)
               endif
               Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(4)%Obs_Latt(imj,1,no_I,no_J) &
               & + ZP*ZS * spin2(GR, GRC, ZN, I1, J1)
               Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) = Obs_eq(5)%Obs_Latt(imj,1,no_I,no_J) &
               & + ZP*ZS * spin3(GR, GRC, ZN, I1, J1)

               ! Dimer  [dimer(GR, GRC, NSUN, i, j, k, l)]
              !  if ( no_I /= 1 .or. no_J /= 1 ) then
              !     write(error_unit,*) "Error in Dimer obeservation: Not expected orbital", no_I, no_J
              !     error stop
              !  endif
               do n1=1,2  ! n1=1: Horizontal dimer, n1=2: Vertical dimer
                nI = n1 + (no_I-1)*2
               do no2_I2=1, N_Spin
                  if ( n1 == 1 .and. latt_id < 3 ) I2 = Invlist2( latt%nnlist(I,1,0), no_I, no2_I2 )
                  if ( n1 == 2 .and. latt_id == 1 ) I2 = Invlist2( latt%nnlist(I,0,1), no_I, no2_I2 )
                  if ( n1 == 2 .and. latt_id == 2 ) then
                     if (no_I < latt_unit%norb) then
                        I2 = Invlist2( I, no_I+1, no2_I2 )
                     else
                        I2 = Invlist2( I, 1, no2_I2 )
                     endif
                  endif
                  if ( latt_id == 3 ) then  ! Open_square
                     I1x = list_orb(I, 1)
                     I1y = list_orb(I, 2)
                     if ( n1 == 1 ) then
                        I2x = 1
                        if (I1x < L1) I2x = I1x+1
                        I2y = I1y
                     else if ( n1 == 2 ) then
                        I2y = 1
                        if (I1y < L1 ) I2y = I1y+1
                        I2x = I1x
                     endif
                     no_I2 = Invlist_orb(I2x, I2y)
                     I2 = Invlist2( I, no_I2, no2_I2 )
                  endif
                  do n2=1,2  ! n2=1: Horizontal dimer, n2=2: Vertical dimer
                    nJ = n2 + (no_J-1)*2
                  do no2_J2=1, N_Spin
                     if ( n2 == 1 .and. latt_id < 3 ) J2 = Invlist2( latt%nnlist(J,1,0), no_J, no2_J2 )
                     if ( n2 == 2 .and. latt_id == 1 ) J2 = Invlist2( latt%nnlist(J,0,1), no_J, no2_J2 )
                     if ( n2 == 2 .and. latt_id == 2 ) then
                        if (no_J < latt_unit%norb) then
                           J2 = Invlist2( J, no_J+1, no2_J2 )
                        else
                           J2 = Invlist2( J, 1, no2_J2 )
                        endif
                     endif
                     if ( latt_id == 3 ) then  ! Open_square
                        J1x = list_orb(J, 1)
                        J1y = list_orb(J, 2)
                        if ( n2 == 1 ) then
                           J2x = 1
                           if (J1x < L1) J2x = J1x+1
                           J2y = J1y
                        else if ( n2 == 2 ) then
                           J2y = 1
                           if (J1y < L1 ) J2y = J1y+1
                           J2x = J1x
                        endif
                        no_J2 = Invlist_orb(J2x, J2y)
                        J2 = Invlist2( I, no_J2, no2_J2 )
                     endif
                     Obs_eq(6)%Obs_Latt(imj,1,nI,nJ) = Obs_eq(6)%Obs_Latt(imj,1,nI,nJ) &
                     & + ZP*ZS * DDDD(GR, GRC, ZN, I1, I2, J1, J2) *cmplx( 0.25d0,  0.d0, kind(0.D0) )
                     Obs_eq(7)%Obs_Latt(imj,1,nI,nJ) = Obs_eq(7)%Obs_Latt(imj,1,nI,nJ) &
                     & + ZP*ZS * dimer2(GR, GRC, ZN, I1, I2, J1, J2)
                     Obs_eq(8)%Obs_Latt(imj,1,nI,nJ) = Obs_eq(8)%Obs_Latt(imj,1,nI,nJ) &
                     & + ZP*ZS * dimer3(GR, GRC, ZN, I1, I2, J1, J2)
                     Obs_eq(9)%Obs_Latt(imj,1,nI,nJ) = Obs_eq(9)%Obs_Latt(imj,1,nI,nJ) &
                     & + ZP*ZS * dimer_z(GR, GRC, I1, I2, J1, J2)
                  enddo
                  enddo
               enddo
               enddo
            ENDDO
            ! Den0
            Obs_eq(2)%Obs_Latt0(no_I) = Obs_eq(2)%Obs_Latt0(no_I) + ZN * GRC(I1,I1,1) * ZP * ZS

            ! Dimer0
            do n1=1,2
              nI = n1 + (no_I-1)*2
            do no2_I2=1, N_Spin
               if ( n1 == 1 .and. latt_id < 3  ) I2 = Invlist2( latt%nnlist(I,1,0), no_I, no2_I2 )
               if ( n1 == 2 .and. latt_id == 1 ) I2 = Invlist2( latt%nnlist(I,0,1), no_I, no2_I2 )
               if ( n1 == 2 .and. latt_id == 2 ) then
                  if (no_I < latt_unit%norb) then
                     I2 = Invlist2( I, no_I+1, no2_I2 )
                  else
                     I2 = Invlist2( I, 1, no2_I2 )
                  endif
               endif
               if ( latt_id == 3 ) then  ! Open_square
                  I1x = list_orb(I, 1)
                  I1y = list_orb(I, 2)
                  if ( n1 == 1 ) then
                     I2x = 1
                     if (I1x < L1) I2x = I1x+1
                     I2y = I1y
                  else if ( n1 == 2 ) then
                     I2y = 1
                     if (I1y < L1 ) I2y = I1y+1
                     I2x = I1x
                  endif
                  no_I2 = Invlist_orb(I2x, I2y)
                  I2 = Invlist2( I, no_I2, no2_I2 )
               endif
               Obs_eq(6)%Obs_Latt0(nI) = Obs_eq(6)%Obs_Latt0(nI) &
                  & + ZP*ZS * DD(GR, GRC, ZN, I1, I2) * cmplx( 0.5d0,  0.d0, kind(0.D0) )
               Obs_eq(7)%Obs_Latt0(nI) = Obs_eq(7)%Obs_Latt0(nI) &
                  & + ZP*ZS * spin2(GR, GRC, ZN, I1, I2)
               Obs_eq(8)%Obs_Latt0(nI) = Obs_eq(8)%Obs_Latt0(nI) &
                  & + ZP*ZS * dimer3_0(GR, GRC, ZN, I1, I2)
               Obs_eq(9)%Obs_Latt0(nI) = Obs_eq(9)%Obs_Latt0(nI) &
                  & + ZP*ZS * spin_z(GR, GRC, I1, I2)
               ! bonds
               Obs_eq(10)%Obs_Latt(I,1,nI,nI) = Obs_eq(10)%Obs_Latt(I,1,nI,nI) &
                  & + ZP*ZS * spin3(GR, GRC, ZN, I1, I2) * latt%N
               Obs_eq(11)%Obs_Latt(I,1,nI,nI) = Obs_eq(11)%Obs_Latt(I,1,nI,nI) &
                  & + ZP*ZS * dimer3_0(GR, GRC, ZN, I1, I2) * latt%N
               ! bonds_back
               Z_back(nI) = Z_back(nI) + spin3(GR, GRC, ZN, I1, I2)
            enddo
            enddo
            
         ENDDO
         Z_back = Z_back / latt%N
         call Obs_scal(4)%measure(Z_back, Phase )

      end Subroutine Obser_AFM_N_S


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
      Subroutine ObserT_AFM_N_S(NT, GT0,G0T,G00,GTT, PHASE, Mc_step_weight)
         Implicit none

         Integer         , INTENT(IN) :: NT
         Complex (Kind=Kind(0.d0)), INTENT(IN) :: GT0(Ndim,Ndim,N_FL),G0T(Ndim,Ndim,N_FL)
         Complex (Kind=Kind(0.d0)), INTENT(IN) :: G00(Ndim,Ndim,N_FL),GTT(Ndim,Ndim,N_FL)
         Complex (Kind=Kind(0.d0)), INTENT(IN) :: Phase
         Real    (Kind=Kind(0.d0)), INTENT(IN) :: Mc_step_weight

         !Locals
         Complex (Kind=Kind(0.d0)) :: ZN, ZP, ZS, ZS0_sq
         Integer :: IMJ, I, J, I1, J1, no_I, no_J

         ZS0_sq = cmplx( ( dble(N_SUN)**2 + dble(N_SUN) )/8.d0 , 0.d0, kind(0.D0))

         ZP = PHASE/Real(Phase, kind(0.D0))
         ZS = Real(Phase, kind(0.D0))/Abs(Real(Phase, kind(0.D0)))
         If (NT == 0 ) then
            DO I = 1,Size(Obs_tau,1)
               Obs_tau(I)%N = Obs_tau(I)%N + 1
               Obs_tau(I)%Ave_sign = Obs_tau(I)%Ave_sign + Real(ZS,kind(0.d0))
            ENDDO
         endif

         ZN = cmplx(dble(N_SUN),0.d0, kind(0.D0))
         Do I1 = 1,Ndim
            I    = List2(I1,1)
            no_I = List2(I1,2)
            Do J1 = 1,Ndim
               J    = List2(J1,1)
               no_J = List2(J1,2)
               imj = latt%imj(I,J)
               ! Green
               Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(1)%Obs_Latt(imj,nt+1,no_I,no_J) &
                  & + ZN * GT0(I1,J1,1) * ZP* ZS

               ! Den
               Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(2)%Obs_Latt(imj,nt+1,no_I,no_J) &
                  & + ( ZN*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1))  &
                  &       *(cmplx(1.d0,0.d0,kind(0.d0)) - G00(J1,J1,1))  &
                  &     -  GT0(I1,J1,1)*G0T(J1,I1,1)                       ) * ZN * ZP*ZS

               ! Spin
               if ( I1 == J1 .and. nt == 0 ) then
                  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) &
                  & + ZP*ZS * ZS0_sq
               else
                  Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) = Obs_tau(3)%Obs_Latt(imj,nt+1,no_I,no_J) &
                  & + ZP*ZS * spin_tau(GT0, G0T, G00, GTT, ZN, I1, J1)
               endif
            Enddo

            ! Den0
            Obs_tau(2)%Obs_Latt0(no_I) = Obs_tau(2)%Obs_Latt0(no_I) + &
               &         ZN*(cmplx(1.d0,0.d0,kind(0.d0)) - GTT(I1,I1,1)) * ZP * ZS
         Enddo

      end Subroutine ObserT_AFM_N_S


! !--------------------------------------------------------------------
! !> @author
! !> ALF Collaboration
! !>
! !> @brief
! !> Global moves
! !>
! !> @details
! !>  This routine generates a
! !>  global update  and returns the propability T0_Proposal_ratio  =  T0( sigma_out-> sigma_in ) /  T0( sigma_in -> sigma_out)
! !> @param [IN] nsigma_old,  Type(Fields)
! !> \verbatim
! !>  Old configuration. The new configuration is stored in nsigma.
! !> \endverbatim
! !> @param [OUT]  T0_Proposal_ratio Real
! !> \verbatimam
! !>  T0_Proposal_ratio  =  T0( sigma_new -> sigma_old ) /  T0( sigma_old -> sigma_new)
! !> \endverbatim
! !> @param [OUT]  Size_clust Real
! !> \verbatim
! !>  Size of cluster that will be flipped.
! !> \endverbatim
! !-------------------------------------------------------------------
!       ! Functions for Global moves.  These move are not implemented in this example.
!       Subroutine Global_move_AFM_N_S(,T0_Proposal_ratio,nsigma_old,size_clust)
!          use Random_Wrap

!          Implicit none
!          Real (Kind=Kind(0.d0)), intent(out) :: T0_Proposal_ratio, size_clust
!          Type (Fields),  Intent(IN)  :: nsigma_old

!          Integer :: n_F1, n_F2, n, nt, I1, I2, I3, I4, nc
!          Integer :: offset, last
!          Logical :: bool

!          T0_Proposal_ratio = 1.d0

!          if ( RANF_WRAP() > 0.5d0 ) then
!             bool = .true.
!             size_clust = 1.d0
!          else
!             bool = .false.
!             size_clust = -1.d0
!          endif

!          offset = Ndim + 2*B_intra*L_intra
!          last   = Ndim + 2*B_intra*L_intra + 4*N_Fam*L_FAM

!          nc=0
!          print*, "offset", offset
!          print*, "last", last

!          do n_F1 = 0,N_Fam-1
!             if ( n_F1 == 0 ) then
!                n_F2 = N_Fam-1
!             else
!                n_F2 = n_F1-1
!             endif
!             do n = 1,2*L_FAM
!                I1 = offset + n_F1 * 2*L_FAM + n
!                I2 = offset + n_F2 * 2*L_FAM + n
!                I3 = last - n_F1 * 2*L_FAM - n + 1
!                I4 = last - n_F2 * 2*L_FAM - n + 1
!                print*, I1, I2, Op_V(I1,1)%I_Fam, Op_V(I2,1)%I_Fam
!                print*, I3, I4, Op_V(I3,1)%I_Fam, Op_V(I4,1)%I_Fam
!                do nt = 1,Ltrot
!                   if ( bool ) then
!                      nsigma%f(I1,nt) = nsigma_old%f(I3,nt)
!                      nsigma%f(I2,nt) = nsigma_old%f(I4,nt)
!                   else
!                      nsigma%f(I3,nt) = nsigma_old%f(I1,nt)
!                      nsigma%f(I4,nt) = nsigma_old%f(I2,nt)
!                   endif
!                enddo
!                nc = nc+1
!             enddo
!          enddo
!          print*, nc

!       End Subroutine Global_move_AFM_N_S


      Subroutine get_param_AFM_N_S(attr_double, attr_names_double, &
                           attr_int, attr_names_int, &
                           attr_str, attr_names_str, &
                           attr_logical, attr_names_logical)
         Implicit none
         Real (Kind=Kind(0.d0)), allocatable, intent(out) :: attr_double(:)
         Integer,                allocatable, intent(out) :: attr_int(:)
         Character (len=64),     allocatable, intent(out) :: attr_str(:)
         Logical,                allocatable, intent(out) :: attr_logical(:)
         Character (len=64),     allocatable, intent(out) :: attr_names_double(:), attr_names_int(:)
         Character (len=64),     allocatable, intent(out) :: attr_names_str(:), attr_names_logical(:)

         allocate( attr_int(6), attr_names_int(6) )
         allocate( attr_str(2), attr_names_str(2) )
         allocate( attr_logical(2), attr_names_logical(2) )
         allocate( attr_double(7), attr_names_double(7) )
         attr_int(1) = L1; attr_names_int(1) = 'L1'
         attr_int(2) = L2; attr_names_int(2) = 'L2'
         attr_int(3) = N_SUN; attr_names_int(3) = 'N_SUN'
         attr_int(4) = N_Spin; attr_names_int(4) = 'N_Spin'
         attr_int(5) = pinning_x; attr_names_int(5) = 'pinning_x'
         attr_int(6) = pinning_y; attr_names_int(6) = 'pinning_y'
         attr_str(1) = Lattice_type; attr_names_str(1) ='Lattice_type'
         attr_str(2) = 'AFM_N_S'; attr_names_str(2) ='model_name'
         attr_logical(1) = Symm; attr_names_logical(1) = 'Symm'
         attr_logical(2) = Projector; attr_names_logical(2) = 'Projector'
         attr_double(1) = ham_J; attr_names_double(1) = 'ham_J'
         attr_double(2) = ham_U; attr_names_double(2) = 'ham_U'
         attr_double(3) = ham_J2; attr_names_double(3) = 'ham_J2'
         attr_double(4) = dtau; attr_names_double(4) = 'dtau'
         attr_double(5) = beta; attr_names_double(5) = 'beta'
         attr_double(6) = theta; attr_names_double(6) = 'theta'
         attr_double(7) = pinning_factor; attr_names_double(7) = 'pinning_factor'
      end Subroutine get_param_AFM_N_S


!####################################################################

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

   end submodule Ham_AFM_N_S_mod
