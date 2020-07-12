module compute_flux_module

  use amrex_fort_module, only : amrex_real
  implicit none

  private

  public :: compute_con_flux, compute_diff_flux, compute_ad_flux, compute_source_flux

contains

!----------------------------------------------------------------------------------------
! Subroutine to compute the convective fluxes
subroutine compute_con_flux(level, ddir, nc, lo, hi,  & 
          &                 u0, u0_lo, u0_hi,       &
          &                 flxx, fx_lo, fx_hi,       &
          &                 flxy, fy_lo, fy_hi  )

    integer, intent(in) :: level, ddir, nc, lo(3), hi(3)
    integer, intent(in) :: u0_lo(3), u0_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    real(amrex_real), intent(in   ) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:nc-1)
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:nc-1)

    integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
    real(amrex_real), parameter :: half = 0.5_amrex_real
    integer :: i, j, k, n, ng, nghost(2)
    integer :: ilo, ihi, jlo, jhi, klo, khi


    flxx(:,:,:,pre) = 0.0_amrex_real
    flxy(:,:,:,pre) = 0.0_amrex_real


    ! do k = fx_lo(3), fx_hi(3)
    !   do j = fx_lo(2), fx_hi(2)
    !     do i = fx_lo(1), fx_hi(1)
          
    !     enddo
    !   enddo
    ! enddo
      ! if (ddir == 2) then
      !   flxx = 0.0_amrex_real
      ! else
      do k = fx_lo(3), fx_hi(3)
        do j = fx_lo(2), fx_hi(2)
          do i = fx_lo(1) + 1, fx_hi(1) - 1
            flxx(i,j,k,ro) = half*(u0(i-1,j,k,rou) + u0(i,j,k,rou))
            flxx(i,j,k,rou) = half*( ( (u0(i-1,j,k,rou)**2)/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)**2)/u0(i,j,k,ro) )   )
            flxx(i,j,k,rov) = half*( ( (u0(i-1,j,k,rou)*u0(i-1,j,k,rov))/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )            
            flxx(i,j,k,roE) = half*( ( (u0(i-1,j,k,roE)*u0(i-1,j,k,rou))/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,roE)*u0(i,j,k,rou))/u0(i,j,k,ro) )   )
            if(isnan(flxx(i,j,k,ro)) .or. isnan(flxx(i,j,k,rou)) .or. &
            &  isnan(flxx(i,j,k,rov)) .or. isnan(flxx(i,j,k,roE))) then

                print*,"i = ", i,"j= ",j, " Exiting..NaN found in flxx: ", flxx(i,j,k,ro), ", ", flxx(i,j,k,rou), &
                &      ", ", flxx(i,j,k,rov), ", ", flxx(i,j,k,roE), ", ", flxx(i,j,k,pre) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo

      ! zero-order extrapolation for end points of fx
      do n = 0,nc-1
        do k = fx_lo(3), fx_hi(3)
          do j = fx_lo(2), fx_hi(2)
            flxx(fx_lo(1),j,k,n)  = flxx(fx_lo(1)+1,j,k,n)           
            flxx(fx_hi(1),j,k,n)  = flxx(fx_hi(1)-1,j,k,n)                    
          enddo
        enddo
      enddo
    ! endif

    ! if(ddir == 1) then
    !   flxy = 0.0_amrex_real
    ! else
      do k = fy_lo(3), fy_hi(3)
        do j = fy_lo(2) + 1, fy_hi(2) - 1
          do i = fy_lo(1), fy_hi(1)
            flxy(i,j,k,ro)  = half*( u0(i,j-1,k,rov) + u0(i,j,k,rov) )
            flxy(i,j,k,rou) = half*( ( (u0(i,j-1,k,rou)*u0(i,j-1,k,rov))/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )  
            flxy(i,j,k,rov) = half*( ( (u0(i,j-1,k,rov)**2)/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,rov)**2)/u0(i,j,k,ro) )   )            
            flxy(i,j,k,roE) = half*( ( (u0(i,j-1,k,roE)*u0(i,j-1,k,rov))/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,roE)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )  

            ! print*,"i= ",i,"j= ",j,"flxy= ",flxy(i,j,k,ro), flxy(i,j,k,rou), flxy(i,j,k,rov), flxy(i,j,k,roE)
            if(isnan(flxy(i,j,k,ro)) .or. isnan(flxy(i,j,k,rou)) .or. &
            &  isnan(flxy(i,j,k,rov)) .or. isnan(flxy(i,j,k,roE))) then
                print*,"i = ", i, " Exiting..NaN found in flxy: ", &
                &            flxy(i,j,k,ro), ", ", flxy(i,j,k,rou), &
                &      ", ", flxy(i,j,k,rov), ", ", flxy(i,j,k,roE), ", ", flxy(i,j,k,pre) 
                call exit(123)
            endif           
          enddo
        enddo
      enddo

      ! zero-order extrapolation for end points of fx
      do n = 0,nc-1
        do k = fy_lo(3), fy_hi(3)
          do i = fy_lo(1), fy_hi(1)
            flxy(i,fy_lo(2),k,n) = flxy(i,fy_lo(2)+1,k,n)
            flxy(i,fy_hi(2),k,n) = flxy(i,fy_hi(2)-1,k,n)                 
          enddo
        enddo
      enddo
    ! endif

      print*,"max(abs(flxy))= ", maxval(abs(flxy(:,:,:,rov)))

end subroutine compute_con_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute diffusion stage fluxes 
subroutine compute_diff_flux( level, ddir, nc, dtdx, dtdy, lo, hi,  & 
            &                 u0, u0_lo, u0_hi,      &
            &                 vx, vx_lo, vx_hi,        & 
            &                 vy, vy_lo, vy_hi,         &
            &                 fldx, fdx_lo, fdx_hi,    &
            &                 fldy, fdy_lo, fdy_hi     )

  integer, intent(in) :: level, nc, ddir
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fdx_lo(3), fdx_hi(3)
  integer, intent(in) :: fdy_lo(3), fdy_hi(3)
  real(amrex_real), intent(inout) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: fldx(fdx_lo(1):fdx_hi(1),fdx_lo(2):fdx_hi(2),  &
    &                                     fdx_lo(3):fdx_hi(3),0:nc-2)
  real(amrex_real), intent(inout) :: fldy(fdy_lo(1):fdy_hi(1),fdy_lo(2):fdy_hi(2),  &
    &                                     fdy_lo(3):fdy_hi(3),0:nc-2)


  integer :: i, j, k, n
  real(amrex_real) :: dxdt, dydt, epsx, epsy, nux, nuy
  real(amrex_real), parameter :: one3 = 1.d0/3.d0, one6 = 1.d0/6.d0

  dxdt = 1.d0/dtdx
  dydt = 1.d0/dtdy

  ! calculate the diffusion fluxes in x-direction (do not include the last cell)
  ! if(ddir == 2) then
  !   fldy = 0.0_amrex_real
  ! else
  do k = fdx_lo(3), fdx_hi(3)
    do j = fdx_lo(2), fdx_hi(2)
      do i = fdx_lo(1)+1, fdx_hi(1)-1
        epsx = dtdx*vx(i,j)  
        nux  = one6 + one3*(epsx**2)
        do n = 0,nc-2
          fldx(i,j,k,n) = nux*dxdt*(u0(i-1,j,k,n) - u0(i,j,k,n))    
        enddo
      enddo
    enddo
  enddo

  ! zero-order extrapolation for end points of fldx
  do n = 0,nc-2
    do k = fdx_lo(3), fdx_hi(3)
      do j = fdx_lo(2), fdx_hi(2)
        fldx(fdx_lo(1),j,k,n)  = fldx(fdx_lo(1)+1,j,k,n)            
        fldx(fdx_hi(1),j,k,n)  = fldx(fdx_hi(1)-1,j,k,n)
      enddo                        
    enddo
  enddo
  ! endif

  ! calculate the diffusion fluxes in y-direction (do not include the end cells)
  ! if(ddir == 1) then
  !   fldy = 0.0_amrex_real
  ! else
  do n = 0,nc-2
    do k = fdy_lo(3), fdy_hi(3)
      do j = fdy_lo(2)+1, fdy_hi(2)-1
        do i = fdy_lo(1), fdy_hi(1)
          epsy = dtdy*vy(i,j)
          nuy  = one6 + one3*(epsy**2)
          fldy(i,j,k,n) = nuy*dydt*(u0(i,j-1,k,n) - u0(i,j,k,n))
        enddo 
      enddo
    enddo
  enddo

  ! zero-order extrapolation for end points of fldy
  do n = 0, nc-2
    do k = fdy_lo(3), fdy_hi(3)
      do i = fdy_lo(1), fdy_hi(1)
        fldy(i,fdy_lo(2),k,n) = fldy(i,fdy_lo(2)+1,k,n)
        fldy(i,fdy_hi(2),k,n) = fldy(i,fdy_hi(2)-1,k,n)
      enddo
    enddo
  enddo
! endif

end subroutine compute_diff_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute source term fluxes 
subroutine compute_source_flux( level, nc, ddir, dtdx, dtdy, lo, hi,  & 
            &                 u0, u0_lo, u0_hi,      &
            &                 vx, vx_lo, vx_hi,        & 
            &                 vy, vy_lo, vy_hi,        &
            &                 flsx, fsx_lo, fsx_hi,    &
            &                 flsy, fsy_lo, fsy_hi     )  
  integer, intent(in) :: level, nc, ddir
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fsx_lo(3), fsx_hi(3)
  integer, intent(in) :: fsy_lo(3), fsy_hi(3)
  real(amrex_real), intent(inout) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: flsx(fsx_lo(1):fsx_hi(1),fsx_lo(2):fsx_hi(2),  &
    &                                     fsx_lo(3):fsx_hi(3),0:nc-2)
  real(amrex_real), intent(inout) :: flsy(fsy_lo(1):fsy_hi(1),fsy_lo(2):fsy_hi(2),  &
    &                                     fsy_lo(3):fsy_hi(3),0:nc-2)

  integer :: i, j, k, n
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  ! update source terms in x-direction (we store dP_x, d(Pu), d(Pv))
  flsx(:,:,:,ro) = 0.0_amrex_real
  flsy(:,:,:,ro) = 0.0_amrex_real
  
  if(ddir == 2) then
    flsx = 0.0_amrex_real
  else
  do k = fsx_lo(3), fsx_hi(3)
    do j = fsx_lo(2), fsx_hi(2)
      do i = fsx_lo(1)+1, fsx_hi(1)-1
        flsx(i,j,k,rov) = 0.0_amrex_real
        flsx(i,j,k,rou) = 0.5_amrex_real*( u0(i-1,j,k,pre) + u0(i,j,k,pre) )
        flsx(i,j,k,roE) = 0.5_amrex_real*( (u0(i-1,j,k,pre)*u0(i-1,j,k,rou)/u0(i-1,j,k,ro))  &
        &               + ( u0(i,j,k,pre)*u0(i,j,k,rou)/u0(i,j,k,ro) ) )
      enddo
    enddo
  enddo

  ! zero-order extrapolation for end points of flsx
  do n = rou,roE
    do k = fsx_lo(3), fsx_hi(3)
      do j = fsx_lo(2), fsx_hi(2)
        flsx(fsx_lo(1),j,k,n)  = flsx(fsx_lo(1)+1,j,k,n)            
        flsx(fsx_hi(1),j,k,n)  = flsx(fsx_hi(1)-1,j,k,n)
      enddo                        
    enddo
  enddo
endif

  ! update source terms corresponding to y-direction
  ! if(ddir == 1) then
  !   flsy = 0.0_amrex_real
  ! else
  do k = fsy_lo(3), fsy_hi(3)
    do j = fsy_lo(2)+1,fsy_hi(2)-1
      do i = fsy_lo(1),fsy_hi(1)
        flsy(i,j,k,rou) = 0.0_amrex_real
        flsy(i,j,k,rov) = 0.5_amrex_real*(u0(i,j-1,k,pre) + u0(i,j,k,pre))
        flsy(i,j,k,roE) = 0.5_amrex_real*( (u0(i,j-1,k,pre)*u0(i,j-1,k,rov)/u0(i,j-1,k,ro)) &
        &               + ( u0(i,j,k,pre)*u0(i,j,k,rov)/u0(i,j,k,ro) ) )
        ! print*,"i= ",i,"j= ",j,"flsy= ",flsy(i,j,k,ro), flsy(i,j,k,rou), flsy(i,j,k,rov), flsy(i,j,k,roE)
      enddo
    enddo
  enddo

  ! zero order extrapolation to end points of flsy
  do n = rou, roE 
    do k = fsy_lo(3),fsy_hi(3)
      do i = fsy_lo(1),fsy_hi(1)
        flsy(i,fsy_lo(2),k,n) = flsy(i,fsy_lo(2)+1,k,n)
        flsy(i,fsy_hi(2),k,n) = flsy(i,fsy_hi(2)-1,k,n)
      enddo
    enddo
  enddo
! endif


end subroutine compute_source_flux
!------------------------------------------------------------------------------------------
! Subroutine to compute anti-diffusive fluxes 
subroutine compute_ad_flux( level, ddir, time, nc, dtdx, dtdy, lo, hi,  & 
      &                 utemp, ut_lo, ut_hi,     &
      &                 ucx, ucx_lo, ucx_hi,     &
      &                 ucy, ucy_lo, ucy_hi,     &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi     )

  integer, intent(in) :: level, nc, ddir
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy, time
  integer, intent(in) :: ut_lo(3), ut_hi(3)
  integer, intent(in) :: ucx_lo(3), ucx_hi(3)
  integer, intent(in) :: ucy_lo(3), ucy_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: ftx_lo(3), ftx_hi(3)
  integer, intent(in) :: fty_lo(3), fty_hi(3)
  real(amrex_real), intent(inout) :: utemp(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
    &                                     ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
    &                                     ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),0:nc-1)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: fltx(ftx_lo(1):ftx_hi(1),ftx_lo(2):ftx_hi(2),  &
    &                                     ftx_lo(3):ftx_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: flty(fty_lo(1):fty_hi(1),fty_lo(2):fty_hi(2),  &
    &                                     fty_lo(3):fty_hi(3),0:nc-1)

  
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  real(amrex_real), parameter :: half = 0.5_amrex_real, one4 = 0.25_amrex_real, &
  &                              one6 = 1.0_amrex_real/6.0_amrex_real
  integer :: i, j, k, n
  integer :: ilo, ihi, jlo, jhi, klo, khi
  real(amrex_real) :: nuxx, nuxy, nuyx, nuyy, epsx, epsy
  ! Following variables are used to compute the conserved variables at the corners which is required 
  ! for the prelimiting step
  real(amrex_real) :: ut1(0:nc-2), ut2(0:nc-2)
  real(amrex_real) :: sgn(0:nc-2), du(0:nc-2), d2u(0:nc-2), fltmp(0:nc-2)
  character(len=128) :: dirchar, fname

  write(dirchar,fmt='(i2.2)') ddir
  ! first calculate x-part of anti-diffusive fluxes (these are at faces normal to x eg. (i+1/2,j) )
  do n = 0,nc-2
    do k = ftx_lo(3), ftx_hi(3)
      do j = ftx_lo(2)+1, ftx_hi(2)-1
        do i = ftx_lo(1)+1, ftx_hi(1)-1
          epsx = dtdx*vx(i,j)  
          epsy = one4*dtdy*(vy(i,j) + vy(i,j+1) + vy(i-1,j) + vy(i-1,j+1))
          nuxx = one6 - one6*(epsx**2)
          nuxy = -half*epsx*epsy
          ! ut1 and ut2 store the values at the corners of the face at which f_ad is calculated
          ut1(n) = one4*(ucx(i,j,k,n) + ucx(i,j+1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j+1,k,n))
          ut2(n) = one4*(ucx(i,j,k,n) + ucx(i,j-1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j-1,k,n))

          fltx(i,j,k,n) = nuxx*(ucx(i,j,k,n) - ucx(i-1,j,k,n)) &
          &             + nuxy*(ut1(n) - ut2(n))

          ! if(j == lo(2) .and. n == nc-2) then
          !   print*,"i= ", i, "fltx= ",fltx(i,j,k,ro),", ",fltx(i,j,k,rou), &
          !   &       ", ",fltx(i,j,k,rov),", ",fltx(i,j,k,roE)
          ! endif

        enddo
      enddo
    enddo
    ! extrapolate to the end points (zero-order extrapolation)
    do k = ftx_lo(3), ftx_hi(3)
      do j = ftx_lo(2), ftx_hi(2)
        fltx(ftx_lo(1),j,k,n) = fltx(ftx_lo(1)+1,j,k,n)
        fltx(ftx_hi(1),j,k,n) = fltx(ftx_hi(1)-1,j,k,n)
      enddo
      do i = ftx_lo(1), ftx_hi(1)
        fltx(i,ftx_lo(2),k,n) = fltx(i,ftx_lo(2)+1,k,n)
        fltx(i,ftx_hi(2),k,n) = fltx(i,ftx_hi(2)-1,k,n)
      enddo
    enddo
  enddo
    !-------------------------------------------------------
    ! print*,"name= ",name, "length= ",length
    ! fname = "fltxd" // trim(dirchar) // ".txt"
    ! open(unit=111,file=fname)
    ! write(111,1100) time
    ! 1100 format('Time= ',F10.5) 
    ! write(111,*) "# i density  x-mom y-mom energy"
    ! ! print*,"lo(1)= ", phi_lo(1), "hi(1)= ",phi_hi(1)
    ! do k = lo(3), hi(3)
    !   do j = lo(2), hi(2)
    !     do i = lo(1), hi(1)
    !       if(j == lo(2)) then
    !         WRITE(111,1200) i, ucx(i,j,k,ro), ucx(i,j,k,rou), ucx(i,j,k,rov), ucx(i,j,k,roE)
    !         1200 format(I5,2x,F14.8,2x,F14.8,2x,F14.8,2x,F14.8)
    !         !-----------------------------------------------------            
    !         ! WRITE(111,1200) i, fltx(i,j,k,ro), fltx(i,j,k,rou), fltx(i,j,k,rov), fltx(i,j,k,roE)
    !         ! 1200 format(I5,2x,F14.8,2x,F14.8,2x,F14.8,2x,F14.8)
    !       endif
    !     enddo
    !   enddo
    ! enddo
    ! close(111)
    !-----------------------------------------------------------
! calculate y-part of anti-diffusive fluxes (these are at faces normal to x eg. (i,j+1/2) )
  do n = 0, nc-2
    do k = fty_lo(3),fty_hi(3)
      do j = fty_lo(2)+1,fty_hi(2)-1
        do i = fty_lo(1)+1,fty_hi(1)-1
          epsy = dtdy*vy(i,j)
          epsx = one4*dtdx*(vx(i,j) + vx(i+1,j) + vx(i,j-1) + vx(i+1,j-1))
          nuyy = one6 - one6*(epsy**2)
          nuxy = -half*epsx*epsy
          ut1(n) = one4*(ucy(i,j,k,n) + ucy(i+1,j,k,n) + ucy(i,j-1,k,n) + ucy(i+1,j-1,k,n))
          ut2(n) = one4*(ucy(i,j,k,n) + ucy(i-1,j,k,n) + ucy(i,j-1,k,n) + ucy(i-1,j-1,k,n))  
          
          flty(i,j,k,n) = nuyy*(ucy(i,j,k,n) - ucy(i,j-1,k,n))  &
          &             + nuxy*(ut1(n) - ut2(n))

          ! if(j == lo(2) .and. n == nc-2) then
          !   print*,"i= ", i, "flty= ",flty(i,j,k,ro),", ",flty(i,j,k,rou), &
          !   &       ", ",flty(i,j,k,rov),", ",flty(i,j,k,roE)
          ! endif        
        enddo
      enddo
    enddo
    ! extrapolate to the end points (zero-order extrapolation)
    do k = fty_lo(3), fty_hi(3)
      do j = fty_lo(2), fty_hi(2)
        flty(fty_lo(1),j,k,n) = flty(fty_lo(1)+1,j,k,n)
        flty(fty_hi(1),j,k,n) = flty(fty_hi(1)-1,j,k,n)
      enddo
      do i = fty_lo(1), fty_hi(1)
        flty(i,fty_lo(2),k,n) = flty(i,fty_lo(2)+1,k,n)
        flty(i,fty_hi(2),k,n) = flty(i,fty_hi(2)-1,k,n)
      enddo
    enddo
  enddo
    !-------------------------------------------------------
    ! print*,"name= ",name, "length= ",length
    ! fname = "fltyd" // trim(dirchar) // ".txt"
    ! open(unit=112,file=fname)
    ! write(112,1100) time
    ! write(112,*) "# j density  x-mom y-mom energy"
    ! ! print*,"lo(1)= ", phi_lo(1), "hi(1)= ",phi_hi(1)
    ! do k = lo(3), hi(3)
    !   do j = lo(2), hi(2)
    !     do i = lo(1), hi(1)
    !       if(i == lo(1)) then
    !         WRITE(112,1201) j, ucy(i,j,k,ro), ucy(i,j,k,rou), ucy(i,j,k,rov), ucy(i,j,k,roE)
    !         1201 format(I5,2x,F14.8,2x,F14.8,2x,F14.8,2x,F14.8)
    !         ! print*,"i=",i,"phi=",phi(i,j,k,1),"x=",x
    !         ! WRITE(112,1201) j, flty(i,j,k,ro), flty(i,j,k,rou), flty(i,j,k,rov), flty(i,j,k,roE)
    !         ! 1201 format(I5,2x,F14.8,2x,F14.8,2x,F14.8,2x,F14.8)
    !       endif
    !     enddo
    !   enddo
    ! enddo
    ! close(112)
    !-----------------------------------------------------------
! Perform the prelimiting step for flux correction (for x-direction fluxes)
do n = 0, nc-2
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1)-1, hi(1)+1
        fltmp(n) = abs(fltx(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, utemp(i,j,k,n) - utemp(i-1,j,k,n))
        du(n)    = utemp(i-1,j,k,n) - utemp(i-2,j,k,n)
        d2u(n)   = utemp(i+1,j,k,n) - utemp(i,j,k,n) 
        fltx(i,j,k,n) = sgn(n)*max( 0.0_amrex_real, min(fltmp(n), sgn(n)*du(n), sgn(n)*d2u(n)) )  
      enddo
    enddo
  enddo
enddo

! Perform the prelimiting step for flux correction (for y-direction fluxes)
do n = 0, nc-2
  do k = lo(3), hi(3)
    do j = lo(2)-1, hi(2)+1
      do i = lo(1), hi(1)
        fltmp(n) = abs(flty(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, utemp(i,j,k,n) - utemp(i,j-1,k,n))
        du(n)    = utemp(i,j-1,k,n) - utemp(i,j-2,k,n)
        d2u(n)   = utemp(i,j+1,k,n) - utemp(i,j,k,n) 
        flty(i,j,k,n) = sgn(n)*max( 0.0_amrex_real, min(fltmp(n), sgn(n)*du(n), sgn(n)*d2u(n)) )  
      enddo
    enddo
  enddo
enddo

! This gives us f^a'(i+1/2,j) and f^a'(i,j+1/2) for all faces 



end subroutine compute_ad_flux


end module compute_flux_module
