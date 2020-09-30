module compute_flux_module

  use amrex_fort_module, only : amrex_real
  implicit none
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, ent = 5
  real(amrex_real), parameter :: half = 0.5_amrex_real, one3 = 1.d0/3.d0, one6 = 1.d0/6.d0

  private

  public :: compute_con_flux, compute_diff_flux, compute_ad_flux, compute_source_flux

contains

!----------------------------------------------------------------------------------------
! Subroutine to compute the convective fluxes
subroutine compute_con_flux(level, nc, lo, hi,  & 
          &                 u0, u0_lo, u0_hi,       &
          &                 flxx, fx_lo, fx_hi,       &
          &                 flxy, fy_lo, fy_hi )

    integer, intent(in) :: level, nc, lo(3), hi(3)
    integer, intent(in) :: u0_lo(3), u0_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    real(amrex_real), intent(in   ) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),ro:ent)
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),ro:ent)
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),ro:ent)

    ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
    ! real(amrex_real), parameter :: half = 0.5_amrex_real
    integer :: i, j, k, n, ng, nghost(2)
    integer :: ilo, ihi, jlo, jhi, klo, khi
    character(len=128) :: fname, dirchar, levchar
    integer :: id(3)

    write(levchar,fmt='(i2.2)') level

    flxx(:,:,:,pre:ent) = 0.0_amrex_real
    flxy(:,:,:,pre:ent) = 0.0_amrex_real

      do k = fx_lo(3), fx_hi(3)
        do j = fx_lo(2), fx_hi(2)
          do i = fx_lo(1) + 1, fx_hi(1) - 1
            flxx(i,j,k,ro)  = half*( u0(i-1,j,k,rou) + u0(i,j,k,rou) )
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
      do n = ro,pre
        do k = fx_lo(3), fx_hi(3)
          do j = fx_lo(2), fx_hi(2)
            flxx(fx_lo(1),j,k,n)  = flxx(fx_lo(1)+1,j,k,n)           
            flxx(fx_hi(1),j,k,n)  = flxx(fx_hi(1)-1,j,k,n)                    
          enddo
        enddo
      enddo

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
      do n = ro,ent
        do k = fy_lo(3), fy_hi(3)
          do i = fy_lo(1), fy_hi(1)
            flxy(i,fy_lo(2),k,n) = flxy(i,fy_lo(2)+1,k,n)
            flxy(i,fy_hi(2),k,n) = flxy(i,fy_hi(2)-1,k,n)                 
          enddo
        enddo
      enddo

end subroutine compute_con_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute diffusion stage fluxes 
subroutine compute_diff_flux( level, nc, dtdx, dtdy, lo, hi,  & 
            &                 u0, u0_lo, u0_hi,      &
            &                 vx, vx_lo, vx_hi,        & 
            &                 vy, vy_lo, vy_hi,         &
            &                 fldx, fdx_lo, fdx_hi,    &
            &                 fldy, fdy_lo, fdy_hi     )

  integer, intent(in) :: level, nc
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fdx_lo(3), fdx_hi(3)
  integer, intent(in) :: fdy_lo(3), fdy_hi(3)
  real(amrex_real), intent(inout) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),ro:ent)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: fldx(fdx_lo(1):fdx_hi(1),fdx_lo(2):fdx_hi(2),  &
    &                                     fdx_lo(3):fdx_hi(3),0:nc-2)
  real(amrex_real), intent(inout) :: fldy(fdy_lo(1):fdy_hi(1),fdy_lo(2):fdy_hi(2),  &
    &                                     fdy_lo(3):fdy_hi(3),0:nc-2)


  integer :: i, j, k, n
  real(amrex_real) :: dxdt, dydt, epsx, epsy, nux, nuy
  ! real(amrex_real), parameter :: one3 = 1.d0/3.d0, one6 = 1.d0/6.d0, half = 0.5_amrex_real

  ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  character(len=128) :: fname, dirchar, levchar

  write(levchar,fmt='(i2.2)') level

  dxdt = 1.d0/dtdx
  dydt = 1.d0/dtdy

  ! calculate the diffusion fluxes in x-direction (do not include the last cell)
  do k = fdx_lo(3), fdx_hi(3)
    do j = fdx_lo(2), fdx_hi(2)
      do i = fdx_lo(1)+1, fdx_hi(1)-1
        ! epsx = dtdx*vx(i,j)  
        ! write epsx in terms of velocity at the beginning of the timestep
        ! get vx using the conserved variables
        epsx = dtdx*half*((u0(i-1,j,k,rou)/u0(i-1,j,k,ro)) + (u0(i,j,k,rou)/u0(i,j,k,ro)))
        nux  = one6 + one3*(epsx**2)
        do n = ro,roE
          fldx(i,j,k,n) = nux*dxdt*(u0(i-1,j,k,n) - u0(i,j,k,n))    
        enddo
      enddo
    enddo
  enddo

  ! zero-order extrapolation for end points of fldx
  do n = ro,roE
    do k = fdx_lo(3), fdx_hi(3)
      do j = fdx_lo(2), fdx_hi(2)
        fldx(fdx_lo(1),j,k,n)  = fldx(fdx_lo(1)+1,j,k,n)            
        fldx(fdx_hi(1),j,k,n)  = fldx(fdx_hi(1)-1,j,k,n)
      enddo                        
    enddo
  enddo
  ! endif

  ! calculate the diffusion fluxes in y-direction (do not include the end cells)
  do n = ro,roE
    do k = fdy_lo(3), fdy_hi(3)
      do j = fdy_lo(2)+1, fdy_hi(2)-1
        do i = fdy_lo(1), fdy_hi(1)
          epsy = dtdy*half*( (u0(i,j-1,k,rov)/u0(i,j-1,k,ro))  &
          &    + (u0(i,j,k,rov)/u0(i,j,k,ro)) )
          nuy  = one6 + one3*(epsy**2)
          fldy(i,j,k,n) = nuy*dydt*(u0(i,j-1,k,n) - u0(i,j,k,n))
        enddo 
      enddo
    enddo
  enddo

  ! zero-order extrapolation for end points of fldy
  do n = ro, roE
    do k = fdy_lo(3), fdy_hi(3)
      do i = fdy_lo(1), fdy_hi(1)
        fldy(i,fdy_lo(2),k,n) = fldy(i,fdy_lo(2)+1,k,n)
        fldy(i,fdy_hi(2),k,n) = fldy(i,fdy_hi(2)-1,k,n)
      enddo
    enddo
  enddo

end subroutine compute_diff_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute source term fluxes 
subroutine compute_source_flux( level, nc, dtdx, dtdy, lo, hi,  & 
            &                 u0, u0_lo, u0_hi,      &
            &                 vx, vx_lo, vx_hi,        & 
            &                 vy, vy_lo, vy_hi,        &
            &                 flsx, fsx_lo, fsx_hi,    &
            &                 flsy, fsy_lo, fsy_hi     )  
  integer, intent(in) :: level, nc
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
    &                                     fsx_lo(3):fsx_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: flsy(fsy_lo(1):fsy_hi(1),fsy_lo(2):fsy_hi(2),  &
    &                                     fsy_lo(3):fsy_hi(3),ro:roE)

  integer :: i, j, k, n
  ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  ! update source terms in x-direction (we store dP_x, d(Pu), d(Pv))
  flsx(:,:,:,ro) = 0.0_amrex_real
  flsy(:,:,:,ro) = 0.0_amrex_real
  
  do k = fsx_lo(3), fsx_hi(3)
    do j = fsx_lo(2), fsx_hi(2)
      do i = fsx_lo(1)+1, fsx_hi(1)-1
        flsx(i,j,k,rov) = 0.0_amrex_real
        flsx(i,j,k,rou) = 0.5_amrex_real*( u0(i-1,j,k,pre) + u0(i,j,k,pre) )
        flsx(i,j,k,roE) = 0.5_amrex_real*( (u0(i-1,j,k,pre)*u0(i-1,j,k,rou)/u0(i-1,j,k,ro))  &
        &               + ( u0(i,j,k,pre)*u0(i,j,k,rou)/u0(i,j,k,ro) ) )
        ! get the face velocities based on the vx and vy multifabs
        ! flsx(i,j,k,roE) = 0.5_amrex_real*( u0(i-1,j,k,pre) + u0(i,j,k,pre) )*vx(i,j)
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
! endif

  ! update source terms corresponding to y-direction
  do k = fsy_lo(3), fsy_hi(3)
    do j = fsy_lo(2)+1,fsy_hi(2)-1
      do i = fsy_lo(1),fsy_hi(1)
        flsy(i,j,k,rou) = 0.0_amrex_real
        flsy(i,j,k,rov) = 0.5_amrex_real*(u0(i,j-1,k,pre) + u0(i,j,k,pre))
        flsy(i,j,k,roE) = 0.5_amrex_real*( (u0(i,j-1,k,pre)*u0(i,j-1,k,rov)/u0(i,j-1,k,ro)) &
        &               + ( u0(i,j,k,pre)*u0(i,j,k,rov)/u0(i,j,k,ro) ) )
        ! flsy(i,j,k,roE) = 0.5_amrex_real*( u0(i,j-1,k,pre) + u0(i,j,k,pre) )*vy(i,j) 
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

end subroutine compute_source_flux
!------------------------------------------------------------------------------------------
! Subroutine to compute anti-diffusive fluxes 
subroutine compute_ad_flux( level, time, nc, dtdx, dtdy, lo, hi,  & 
      &                 utemp, ut_lo, ut_hi,     &
      &                 ucx, ucx_lo, ucx_hi,     &
      &                 ucy, ucy_lo, ucy_hi,     &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi     )

  integer, intent(in) :: level, nc
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy, time
  integer, intent(in) :: ut_lo(3), ut_hi(3)
  integer, intent(in) :: ucx_lo(3), ucx_hi(3)
  integer, intent(in) :: ucy_lo(3), ucy_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: ftx_lo(3), ftx_hi(3)
  integer, intent(in) :: fty_lo(3), fty_hi(3)
  real(amrex_real), intent(inout) :: utemp(ut_lo(1):ut_hi(1),ut_lo(2):ut_hi(2),ut_lo(3):ut_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
    &                                     ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
    &                                     ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),ro:pre)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: fltx(ftx_lo(1):ftx_hi(1),ftx_lo(2):ftx_hi(2),  &
    &                                     ftx_lo(3):ftx_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: flty(fty_lo(1):fty_hi(1),fty_lo(2):fty_hi(2),  &
    &                                     fty_lo(3):fty_hi(3),ro:roE)

  
  ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  ! real(amrex_real), parameter :: half = 0.5_amrex_real, one4 = 0.25_amrex_real, &
  ! &                              one6 = 1.0_amrex_real/6.0_amrex_real
  integer :: i, j, k, n
  integer :: ilo, ihi, jlo, jhi, klo, khi
  real(amrex_real) :: nuxx, nuxy, nuyx, nuyy, epsx, epsy
  ! Following variables are used to compute the conserved variables at the corners which is required 
  ! for the prelimiting step
  real(amrex_real) :: ut1(ro:roE), ut2(ro:roE)
  real(amrex_real) :: sgn(ro:roE), du(ro:roE), d2u(ro:roE), fltmp(ro:roE)
  ! real(amrex_real) :: ut1(0:nc-1), ut2(0:nc-1)
  ! real(amrex_real) :: sgn(0:nc-1), du(0:nc-1), d2u(0:nc-1), fltmp(0:nc-1)
  character(len=128) :: dirchar, fname

  ! first calculate x-part of anti-diffusive fluxes (these are at faces normal to x eg. (i+1/2,j) )
  ! print*,"ftx_lo= ",ftx_lo, ", ucx_lo = ",ucx_lo
  ! print*,"ftx_hi= ",ftx_hi, ", ucx_hi= ",ucx_hi
  do n = ro,roE
    do k = ftx_lo(3), ftx_hi(3)
      do j = ftx_lo(2)+1, ftx_hi(2)-1
        do i = ftx_lo(1)+1, ftx_hi(1)-1
          epsx = dtdx*vx(i,j)  
          ! epsy = one4*dtdy*(vy(i,j) + vy(i,j+1) + vy(i-1,j) + vy(i-1,j+1))

          ! epsx = dtdx*half*( (utemp(i-1,j,k,rou)/utemp(i-1,j,k,ro)) &
          ! &    + (utemp(i,j,k,rou)/utemp(i,j,k,ro)) )
          nuxx = one6 - one6*(epsx**2)
          ! nuxy = -half*epsx*epsy
          ! ut1 and ut2 store the values at the corners of the face at which f_ad is calculated
          ! ut1(n) = one4*(ucx(i,j,k,n) + ucx(i,j+1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j+1,k,n))
          ! ut2(n) = one4*(ucx(i,j,k,n) + ucx(i,j-1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j-1,k,n))

          fltx(i,j,k,n) = nuxx*(ucx(i,j,k,n) - ucx(i-1,j,k,n)) 
          ! &
          ! &             + nuxy*(ut1(n) - ut2(n))

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
  
! calculate y-part of anti-diffusive fluxes (these are at faces normal to x eg. (i,j+1/2) )
  do n = ro, roE
  ! do n = 0, nc-1
    do k = fty_lo(3),fty_hi(3)
      do j = fty_lo(2)+1,fty_hi(2)-1
        do i = fty_lo(1)+1,fty_hi(1)-1
          epsy = dtdy*vy(i,j)
          ! epsx = one4*dtdx*(vx(i,j) + vx(i+1,j) + vx(i,j-1) + vx(i+1,j-1))
          nuyy = one6 - one6*(epsy**2)
          ! nuxy = -half*epsx*epsy
          ! ut1(n) = one4*(ucy(i,j,k,n) + ucy(i+1,j,k,n) + ucy(i,j-1,k,n) + ucy(i+1,j-1,k,n))
          ! ut2(n) = one4*(ucy(i,j,k,n) + ucy(i-1,j,k,n) + ucy(i,j-1,k,n) + ucy(i-1,j-1,k,n))  
          
          flty(i,j,k,n) = nuyy*(ucy(i,j,k,n) - ucy(i,j-1,k,n))  
          ! &
          ! &             + nuxy*(ut1(n) - ut2(n))

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
  ! print*,"calculated y-antidiffusive fluxes"
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
if(level == 0) then
  ilo = lo(1)-1; ihi = hi(1)+1
  jlo = lo(2); jhi = hi(2)
  klo = lo(3); khi = hi(3)
else
  ilo = lo(1)-2; ihi = hi(1)+3
  jlo = lo(2)-3; jhi = hi(2)+3
  klo = lo(3); khi = hi(3)
endif
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        ! print*,"i= ",i,"j= ",j,"k= ",k,"n= ",n,"reached here"
        fltmp(n) = abs(fltx(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, utemp(i,j,k,n) - utemp(i-1,j,k,n))
        du(n)    = utemp(i-1,j,k,n) - utemp(i-2,j,k,n)
        d2u(n)   = utemp(i+1,j,k,n) - utemp(i,j,k,n) 
        fltx(i,j,k,n) = sgn(n)*max( 0.0_amrex_real, min(fltmp(n), sgn(n)*du(n), sgn(n)*d2u(n)) )  
      enddo
    enddo
  enddo
  ! extrapolate to the end points (zero-order extrapolation)
  if(level > 0) then
    do k = ftx_lo(3), ftx_hi(3)
      do j = ftx_lo(2), ftx_hi(2)
        fltx(ftx_lo(1)+1,j,k,n) = fltx(ftx_lo(1)+2,j,k,n)
        fltx(ftx_lo(1),j,k,n) = fltx(ftx_lo(1)+2,j,k,n)
        fltx(ftx_hi(1)-1,j,k,n) = fltx(ftx_hi(1)-2,j,k,n)
        fltx(ftx_hi(1),j,k,n) = fltx(ftx_hi(1)-2,j,k,n)
      enddo
      do i = ftx_lo(1), ftx_hi(1)
        fltx(i,ftx_lo(2),k,n) = fltx(i,ftx_lo(2)+1,k,n)
        fltx(i,ftx_hi(2),k,n) = fltx(i,ftx_hi(2)-1,k,n)
      enddo
    enddo
  endif
enddo
! print*,"prelimited x-antidiffusive fluxes"

if(level == 0) then
  ilo = lo(1); ihi = hi(1)
  jlo = lo(2)-1; jhi = hi(2)+1
  klo = lo(3); khi = hi(3)
else
  ilo = lo(1)-3; ihi = hi(1)+3
  jlo = lo(2)-2; jhi = hi(2)+3
  klo = lo(3); khi = hi(3)
endif
! Perform the prelimiting step for flux correction (for y-direction fluxes)
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        fltmp(n) = abs(flty(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, utemp(i,j,k,n) - utemp(i,j-1,k,n))
        du(n)    = utemp(i,j-1,k,n) - utemp(i,j-2,k,n)
        d2u(n)   = utemp(i,j+1,k,n) - utemp(i,j,k,n) 
        flty(i,j,k,n) = sgn(n)*max( 0.0_amrex_real, min(fltmp(n), sgn(n)*du(n), sgn(n)*d2u(n)) )  
      enddo
    enddo
  enddo
! extrapolate to the end points (zero-order extrapolation)
  if(level > 0) then
    do k = fty_lo(3), fty_hi(3)
      do j = fty_lo(2), fty_hi(2)
        flty(fty_lo(1),j,k,n) = flty(fty_lo(1)+1,j,k,n)
        flty(fty_hi(1),j,k,n) = flty(fty_hi(1)-1,j,k,n)
      enddo
      do i = fty_lo(1), fty_hi(1)
        flty(i,fty_lo(2)+1,k,n) = flty(i,fty_lo(2)+2,k,n)
        flty(i,fty_lo(2),k,n) = flty(i,fty_lo(2)+1,k,n)
        flty(i,fty_hi(2),k,n) = flty(i,fty_hi(2)-1,k,n)
      enddo
    enddo
  endif
enddo
! print*,"prelimited y-antidiffusive fluxes"

! This gives us f^a'(i+1/2,j) and f^a'(i,j+1/2) for all faces 



end subroutine compute_ad_flux


end module compute_flux_module
