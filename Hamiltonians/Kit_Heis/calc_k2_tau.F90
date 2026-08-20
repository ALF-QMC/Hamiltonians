!  Copyright (C) 2026 The ALF project
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
!     along with Foobar.  If not, see http://www.gnu.org/licenses/.
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

       Program calc_k2_tau

!--------------------------------------------------------------------
!> @author 
!> ALF-project
!
!> @brief 
!> Calculate time-displaced correlation functions for the term k2 of the magnetotropic susceptibility k.
!
!--------------------------------------------------------------------
         use iso_fortran_env, only: output_unit, error_unit

         Use Errors
         Use MyMats
         use Lattices_v3
         Use Matrix
         use ana_mod

         Implicit none


         Integer :: Nunit, Norb, N_auto, N_BZ_Zones
         logical :: Extended_Zone, projector
         Integer :: nb, no, no1, n, nbins, n_skip, NT, NT1, Lt, N_rebin, N_cov, N_Back
         Integer :: Lt_eff
         Complex (Kind=Kind(0.d0)), allocatable :: Xmean(:), Xcov(:,:)
         Real    (Kind=Kind(0.d0)), parameter :: Zero=1.D-8, pi=acos(-1.d0)
         Real    (Kind=Kind(0.d0)), allocatable :: Phase(:)
         Character (len=64) :: File_out


         Complex (Kind=Kind(0.d0)), allocatable :: SS(:,:,:,:),kmc2t(:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: SD(:,:,:,:),kmc2tph(:,:)
         Real    (Kind=Kind(0.d0)):: Ham_J,Ham_J3,Ham_alphax,Ham_alphay,Ham_alphaz,Ham_Gx
         Real    (Kind=Kind(0.d0)):: Ham_Gy,Ham_Gz,Hab,Htheta,Hphi,Ham_U,dtau,beta
         Real    (Kind=Kind(0.d0)):: Ham_Gx_p,Ham_Gy_p,Ham_Gz_p,Theta
         Real    (Kind=Kind(0.d0)):: n0(1:3),m0(1:3),B1,B2,B3,g1,g2,g3,e1,e2,e3,gA,gC


         Complex (Kind=Kind(0.d0)), allocatable :: bins_raw(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_xx(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_yy(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_zz(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_xy(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_yx(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_xz(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_zx(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_yz(:,:,:,:,:)
         Complex (Kind=Kind(0.d0)), allocatable :: Bins_zy(:,:,:,:,:)
         Type (Lattice)   :: Latt

         NAMELIST /VAR_errors/   n_skip, N_rebin, N_Cov, N_Back, N_auto, N_BZ_Zones, Extended_Zone
         NAMELIST /VAR_Kit_Heis/ Ham_J,Ham_J3,Ham_alphax,Ham_alphay,Ham_alphaz, &
                 Ham_Gx,Ham_Gy,Ham_Gz,Ham_Gx_p,Ham_Gy_p,Ham_Gz_p,Hab,Htheta,Hphi,Ham_U,Dtau,Beta,theta,projector

         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read')
         READ(5,NML=VAR_errors)
         CLOSE(5)
         
         OPEN(UNIT=5,FILE='parameters',STATUS='old',ACTION='read')
         READ(5,NML=VAR_Kit_Heis)
         CLOSE(5)

         g1=2.3d0
         g2=2.3d0
         g3=1.3d0

         !ac plane
         B1=     Hab*sin(Htheta*pi)/(6d0)**0.5d0+Hab*cos(Htheta*pi)/(3d0)**0.5d0
         B2=     Hab*sin(Htheta*pi)/(6d0)**0.5d0+Hab*cos(Htheta*pi)/(3d0)**0.5d0
         B3=-2d0*Hab*sin(Htheta*pi)/(6d0)**0.5d0+Hab*cos(Htheta*pi)/(3d0)**0.5d0
         e1=-1d0/2d0**0.5d0
         e2= 1d0/2d0**0.5d0
         e3= 0d0
         gA=1d0/3d0*(g3+2d0*g1)
         gC=1d0/3d0*(g3-g1)
         n0(1)=-(e2*B3-e3*B2)
         n0(2)=-(e3*B1-e1*B3)
         n0(3)=-(e1*B2-e2*B1)
         m0(1)= e2*(e1*B2-e2*B1) -e3*(e3*B1-e1*B3)
         m0(2)= e3*(e2*B3-e3*B2) -e1*(e1*B2-e2*B1)
         m0(3)= e1*(e3*B1-e1*B3) -e2*(e2*B3-e3*B2)

         call read_obs("Spinxx_tau", bins_raw, Latt, 0)
         Nunit = size(bins_raw, 1)
         LT = size(bins_raw, 2)
         Norb = size(bins_raw, 3)
         Nbins = size(bins_raw, 5)
         Write(output_unit,"(A,I0)") "# of bins:           ", Nbins
         Nbins  = Nbins - n_skip
         Write(output_unit,"(A,I0)") "Effective # of bins: ", Nbins
         Write(output_unit,"(A,I0)") "After rebinnining:   ", Nbins / N_rebin
         if(Nbins / N_rebin < 2) then
           write (error_unit,*) "Effective # of bins smaller than 2. Analysis impossible!"
           stop 1
         endif
         deallocate(bins_raw)

         call read_obs("Spinxx_tau", bins_xx, Latt, n_skip, phase)
         call read_obs("Spinyy_tau", bins_yy, Latt, n_skip)
         call read_obs("Spinzz_tau", bins_zz, Latt, n_skip)
         call read_obs("Spinxy_tau", bins_xy, Latt, n_skip)
         call read_obs("Spinyx_tau", bins_yx, Latt, n_skip)
         call read_obs("Spinxz_tau", bins_xz, Latt, n_skip)
         call read_obs("Spinzx_tau", bins_zx, Latt, n_skip)
         call read_obs("Spinyz_tau", bins_yz, Latt, n_skip)
         call read_obs("Spinzy_tau", bins_zy, Latt, n_skip)

         if (mod(Lt-1,2) == 0 ) then
            Lt_eff = (Lt -1 ) /2   + 1
         else
            Lt_eff = Lt/2 
         endif

         Allocate ( SS(3,3,Lt,Nbins),kmc2t(Lt,Nbins),SD(3,3,Lt,Nbins))
         Allocate (Xmean(Lt_eff), Xcov(Lt_eff,Lt_eff),kmc2tph(Lt_eff,Nbins))

         n = latt%invlistk(0,0)
         SS = cmplx(0.d0,0.d0,kind(0.d0))
         do nb = 1,Nbins
         do no = 1,norb
            do no1 = 1,Norb
               SS(1,1,:,nb) = SS(1,1,:,nb) + bins_xx(n, :, no, no1, nb)*dble(Nunit)/dble(4d0)
               SS(2,2,:,nb) = SS(2,2,:,nb) - bins_yy(n, :, no, no1, nb)*dble(Nunit)/dble(4d0)
               SS(3,3,:,nb) = SS(3,3,:,nb) + bins_zz(n, :, no, no1, nb)*dble(Nunit)/dble(4d0)
               SS(1,2,:,nb) = SS(1,2,:,nb) + 0.5d0*(cmplx(aimag(bins_xy(n, :, no, no1, nb)),-real(bins_xy(n, :, no, no1, nb)),Kind(0.d0)) &
                                                 +cmplx(aimag(bins_yx(n, :, no, no1, nb)),-real(bins_yx(n, :, no, no1, nb)),Kind(0.d0)))*dble(Nunit)/dble(4d0)
               SS(1,3,:,nb) = SS(1,3,:,nb) + 0.5d0*(bins_xz(n, :, no, no1, nb)+bins_zx(n, :, no, no1, nb))*dble(Nunit)/dble(4d0)
               SS(2,3,:,nb) = SS(2,3,:,nb) + 0.5d0*(cmplx(aimag(bins_yz(n, :, no, no1, nb)),-real(bins_yz(n, :, no, no1, nb)),Kind(0.d0)) &
                                                 +cmplx(aimag(bins_zy(n, :, no, no1, nb)),-real(bins_zy(n, :, no, no1, nb)),Kind(0.d0)))*dble(Nunit)/dble(4d0)
            enddo
         enddo
         enddo
         SS(2,1,:,:) = SS(1,2,:,:)
         SS(3,1,:,:) = SS(1,3,:,:)
         SS(3,2,:,:) = SS(2,3,:,:)

         SD = cmplx(0.d0,0.d0,kind(0.d0))
         SD(1,1,:,:)=&
            & gA**2d0*SS(1,1,:,:)&
            &+gC**2d0*(SS(2,2,:,:)+SS(3,3,:,:)+SS(2,3,:,:)+SS(3,2,:,:))&
            &+gA*gC*  (SS(1,2,:,:)+SS(1,3,:,:)+SS(2,1,:,:)+SS(3,1,:,:))

         SD(2,2,:,:)=&
            & gA**2d0*SS(2,2,:,:)&
            &+gC**2d0*(SS(1,1,:,:)+SS(3,3,:,:)+SS(1,3,:,:)+SS(3,1,:,:))&
            &+gA*gC*  (SS(2,1,:,:)+SS(2,3,:,:)+SS(1,2,:,:)+SS(3,2,:,:))

         SD(3,3,:,:)=&
            & gA**2d0*SS(3,3,:,:)&
            &+gC**2d0*(SS(1,1,:,:)+SS(2,2,:,:)+SS(1,2,:,:)+SS(2,1,:,:))&
            &+gA*gC*  (SS(3,1,:,:)+SS(3,2,:,:)+SS(1,3,:,:)+SS(2,3,:,:))

         SD(1,2,:,:)=&
            & gA**2d0*SS(1,2,:,:)&
            &+gC**2d0*(SS(2,1,:,:)+SS(2,3,:,:)+SS(3,1,:,:)+SS(3,3,:,:))&
            &+gA*gC*  (SS(1,1,:,:)+SS(1,3,:,:)+SS(2,2,:,:)+SS(3,2,:,:))

         SD(2,1,:,:)=&
            & gA**2d0*SS(2,1,:,:)&
            &+gC**2d0*(SS(1,2,:,:)+SS(1,3,:,:)+SS(3,2,:,:)+SS(3,3,:,:))&
            &+gA*gC*  (SS(2,2,:,:)+SS(2,3,:,:)+SS(1,1,:,:)+SS(3,1,:,:))

         SD(1,3,:,:)=&
            & gA**2d0*SS(1,3,:,:)&
            &+gC**2d0*(SS(2,1,:,:)+SS(2,2,:,:)+SS(3,1,:,:)+SS(3,2,:,:))&
            &+gA*gC*  (SS(1,1,:,:)+SS(1,2,:,:)+SS(2,3,:,:)+SS(3,3,:,:))

         SD(3,1,:,:)=&
            & gA**2d0*SS(3,1,:,:)&
            &+gC**2d0*(SS(1,2,:,:)+SS(1,3,:,:)+SS(2,2,:,:)+SS(2,3,:,:))&
            &+gA*gC*  (SS(3,2,:,:)+SS(3,3,:,:)+SS(1,1,:,:)+SS(2,1,:,:))

         SD(2,3,:,:)=&
            & gA**2d0*SS(2,3,:,:)&
            &+gC**2d0*(SS(1,1,:,:)+SS(1,2,:,:)+SS(3,1,:,:)+SS(3,2,:,:))&
            &+gA*gC*  (SS(2,1,:,:)+SS(2,2,:,:)+SS(1,3,:,:)+SS(3,3,:,:))

         SD(3,2,:,:)=&
            & gA**2d0*SS(3,2,:,:)&
            &+gC**2d0*(SS(1,1,:,:)+SS(1,3,:,:)+SS(2,1,:,:)+SS(2,3,:,:))&
            &+gA*gC*  (SS(3,1,:,:)+SS(3,3,:,:)+SS(1,2,:,:)+SS(2,2,:,:))

         kmc2t=cmplx(0.d0,0.d0,Kind(0.d0))
         do no = 1,3
            do no1 = 1,3
               kmc2t(:,:)=kmc2t(:,:)+n0(no)*n0(no1)*SD(no,no1,:,:)
            enddo
         enddo

         kmc2tph=cmplx(0.d0,0.d0,Kind(0.d0))
         do nt = 1,Lt_eff
            kmc2tph(nt,:)= kmc2tph(nt,:) + (kmc2t(nt,:)+kmc2t(Lt-nt+1,:))/cmplx(2.d0,0.d0,Kind(0.d0))
         enddo

         call COV(kmc2tph(:,:), phase, Xcov, Xmean, N_rebin )
         File_out = "Spintot"
         Open (Unit=100,File=File_out,status="unknown")
         Write(100,*) Lt_eff,  nbins/N_rebin, real(lt-1,kind(0.d0))*dtau, 1, "PH"
         do nt = 1, Lt_eff
            Write(100,*) &
                  & dble(nt-1)*dtau,  dble(Xmean(nt)), sqrt(abs(dble(Xcov(nt,nt))))
         enddo
         Do nt = 1,LT_eff
            Do nt1 =1,LT_eff
               Write(100,*) dble(Xcov(nt,nt1))
            Enddo
         Enddo
         close(100)
   contains
         subroutine read_obs(name, bins, Latt, n_skip, sgn)
            Character (len=*), intent(in) :: name
            Complex (Kind=Kind(0.d0)), allocatable, intent(out) :: bins(:,:,:,:,:)
            Type (Lattice)                        , intent(out) :: Latt
            Integer, intent(in)  :: n_skip
            Real    (Kind=Kind(0.d0)), allocatable, intent(out), optional :: sgn(:)
            
            Complex (Kind=Kind(0.d0)), pointer :: bins0(:,:), bins_raw(:,:,:,:,:)
            Real    (Kind=Kind(0.d0)), allocatable :: sgn_raw(:)
            Type (Unit_cell)               :: Latt_unit
            Real    (Kind=Kind(0.d0))      :: dtau
            Character (len=:), allocatable :: Channel
            logical :: use_hdf5
            integer :: nb, NBins
#ifdef HDF5
            inquire(file="data.h5", exist=use_hdf5)
            if (use_hdf5) then
               call read_latt_hdf5('data.h5', name, sgn_raw, bins_raw, bins0, Latt, Latt_unit, dtau, Channel)
            else
#endif
               call read_latt(name, sgn_raw, bins_raw, bins0, Latt, Latt_unit, dtau, Channel)
#ifdef HDF5
            endif
#endif
            Nbins = size(bins_raw, 5) - n_skip
            allocate(bins(size(bins_raw,1), size(bins_raw,2), size(bins_raw,3), size(bins_raw,4), Nbins))
            do nb = 1, Nbins
               bins(:,:,:,:,nb) = bins_raw(:,:,:,:,n_skip+nb)
            enddo
            if (present(sgn)) then
               allocate(sgn(Nbins))
               sgn(:) = sgn_raw(n_skip+1:n_skip+Nbins)
            endif
            deallocate(bins0, bins_raw, sgn_raw)
         end subroutine read_obs

       end Program calc_k2_tau
