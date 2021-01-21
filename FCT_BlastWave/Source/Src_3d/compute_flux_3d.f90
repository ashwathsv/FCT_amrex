module compute_flux_module

  use amrex_fort_module, only : amrex_real
  ! use amrex_paralleldescriptor_module
  implicit none
  integer, parameter :: ro = 0, rou = 1, rov = 2, row = 3, roE = 4, pre = 5, mach = 6
  real(amrex_real), parameter :: half = 0.5_amrex_real, one3 = 1.d0/3.d0, one6 = 1.d0/6.d0, one4 = 0.25_amrex_real

  private

  public :: compute_con_flux, compute_diff_flux, compute_ad_flux, prelimit_ad_flux

contains

!----------------------------------------------------------------------------------------
! Subroutine to compute the convective fluxes
subroutine compute_con_flux(level, nc, lo, hi,      & 
          &                 u0, u0_lo, u0_hi,       &
          &                 flxx, fx_lo, fx_hi,     &
          &                 flxy, fy_lo, fy_hi,     &
          &                 flxz, fz_lo, fz_hi      )

    integer, intent(in) :: level, nc, lo(3), hi(3)
    integer, intent(in) :: u0_lo(3), u0_hi(3)
    integer, intent(in) :: fx_lo(3), fx_hi(3)
    integer, intent(in) :: fy_lo(3), fy_hi(3)
    integer, intent(in) :: fz_lo(3), fz_hi(3)
    real(amrex_real), intent(in   ) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),ro:mach)
    real(amrex_real), intent(  out) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),ro:mach)
    real(amrex_real), intent(  out) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),ro:mach)
    real(amrex_real), intent(  out) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),ro:mach)

    integer :: i, j, k, n, ng, nghost(2)
    integer :: ilo, ihi, jlo, jhi, klo, khi
    character(len=128) :: fname, dirchar, levchar
    integer :: id(3)

    write(levchar,fmt='(i2.2)') level

    flxx(:,:,:,pre:mach) = 0.0_amrex_real
    flxy(:,:,:,pre:mach) = 0.0_amrex_real
    flxz(:,:,:,pre:mach) = 0.0_amrex_real

    !$omp parallel do private(k,j,i)
      do k = fx_lo(3), fx_hi(3)
        do j = fx_lo(2), fx_hi(2)
          do i = fx_lo(1) + 1, fx_hi(1) - 1
            flxx(i,j,k,ro)  = half*( u0(i-1,j,k,rou) + u0(i,j,k,rou) )
            flxx(i,j,k,rou) = half*( ( (u0(i-1,j,k,rou)**2)/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)**2)/u0(i,j,k,ro) )   ) &
            &               + half*( u0(i-1,j,k,pre) + u0(i,j,k,pre) )
            flxx(i,j,k,rov) = half*( ( (u0(i-1,j,k,rou)*u0(i-1,j,k,rov))/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )
            flxx(i,j,k,row) = half*( ( (u0(i-1,j,k,rou)*u0(i-1,j,k,row))/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,row))/u0(i,j,k,ro) )   )             
            flxx(i,j,k,roE) = half*( ( (u0(i-1,j,k,roE)*u0(i-1,j,k,rou))/u0(i-1,j,k,ro) ) &
            &               +        ( (u0(i,j,k,roE)*u0(i,j,k,rou))/u0(i,j,k,ro) )   )   &
            &               + half*( ( u0(i-1,j,k,pre)*u0(i-1,j,k,rou)/u0(i-1,j,k,ro) )   &
            &               +        ( u0(i,j,k,pre)*u0(i,j,k,rou)/u0(i,j,k,ro) ) )

            if(isnan(flxx(i,j,k,ro)) .or. isnan(flxx(i,j,k,rou)) .or. &
            &  isnan(flxx(i,j,k,rov)) .or. isnan(flxx(i,j,k,row)) .or. isnan(flxx(i,j,k,roE))) then
                print*,"i = ", i,"j= ",j, " Exiting..NaN found in flxx: ", flxx(i,j,k,ro), ", ", flxx(i,j,k,rou), &
                &      ", ", flxx(i,j,k,rov), ", ", flxx(i,j,k,row), ", ", flxx(i,j,k,roE), ", ", flxx(i,j,k,pre) 
                print*,"i = ", i,"j= ",j, " Exiting..NaN found in u0: ", u0(i,j,k,ro), ", ", u0(i,j,k,rou), &
                &      ", ", u0(i,j,k,rov), ", ", u0(i,j,k,row), ", ", u0(i,j,k,roE), ", ", u0(i,j,k,pre) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo
    !$omp end parallel do

      ! zero-order extrapolation for end points of fy
      do n = ro,pre
        do k = fx_lo(3), fx_hi(3)
          do j = fx_lo(2), fx_hi(2)
            flxx(fx_lo(1),j,k,n)  = flxx(fx_lo(1)+1,j,k,n)           
            flxx(fx_hi(1),j,k,n)  = flxx(fx_hi(1)-1,j,k,n)                    
          enddo
        enddo
      enddo

      !$omp parallel do private(k,j,i)
      do k = fy_lo(3), fy_hi(3)
        do j = fy_lo(2) + 1, fy_hi(2) - 1
          do i = fy_lo(1), fy_hi(1)
            flxy(i,j,k,ro)  = half*( u0(i,j-1,k,rov) + u0(i,j,k,rov) )
            flxy(i,j,k,rou) = half*( ( (u0(i,j-1,k,rou)*u0(i,j-1,k,rov))/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )  
            flxy(i,j,k,rov) = half*( ( (u0(i,j-1,k,rov)**2)/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,rov)**2)/u0(i,j,k,ro) )   ) & 
            &               + half*( u0(i,j-1,k,pre) + u0(i,j,k,pre) )
            flxy(i,j,k,row) = half*( ( (u0(i,j-1,k,rov)*u0(i,j-1,k,row))/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,rov)*u0(i,j,k,row))/u0(i,j,k,ro) )   )            
            flxy(i,j,k,roE) = half*( ( (u0(i,j-1,k,roE)*u0(i,j-1,k,rov))/u0(i,j-1,k,ro) ) &
            &               +        ( (u0(i,j,k,roE)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )   &
            &               + half*( (u0(i,j-1,k,pre)*u0(i,j-1,k,rov)/u0(i,j-1,k,ro)) &
            &               +      ( u0(i,j,k,pre)*u0(i,j,k,rov)/u0(i,j,k,ro) ) )
            if(isnan(flxy(i,j,k,ro)) .or. isnan(flxy(i,j,k,rou)) .or. &
            &  isnan(flxy(i,j,k,rov)) .or. isnan(flxy(i,j,k,row)) .or. isnan(flxy(i,j,k,roE))) then
                print*,"i = ", i, " Exiting..NaN found in flxy: ", &
                &            flxy(i,j,k,ro), ", ", flxy(i,j,k,rou), &
                &      ", ", flxy(i,j,k,rov), ", ", flxy(i,j,k,row), ", ", flxy(i,j,k,roE), ", ", flxy(i,j,k,pre) 
                call exit(123)
            endif           
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! zero-order extrapolation for end points of fy
      do n = ro,mach
        do k = fy_lo(3), fy_hi(3)
          do i = fy_lo(1), fy_hi(1)
            flxy(i,fy_lo(2),k,n) = flxy(i,fy_lo(2)+1,k,n)
            flxy(i,fy_hi(2),k,n) = flxy(i,fy_hi(2)-1,k,n)                 
          enddo
        enddo
      enddo

      !$omp parallel do private(k,j,i)
      do k = fz_lo(3) + 1, fz_hi(3) - 1
        do j = fz_lo(2), fz_hi(2)
          do i = fz_lo(1), fz_hi(1)
            flxz(i,j,k,ro)  = half*( u0(i,j,k-1,row) + u0(i,j,k,row) )
            flxz(i,j,k,rou) = half*( ( (u0(i,j,k-1,row)*u0(i,j,k-1,rou))/u0(i,j,k-1,ro) ) &
            &               +        ( (u0(i,j,k,rou)*u0(i,j,k,row))/u0(i,j,k,ro) )   )  
            flxz(i,j,k,rov) = half*( ( (u0(i,j,k-1,row)*u0(i,j,k-1,rov))/u0(i,j,k-1,ro) ) &
            &               +        ( (u0(i,j,k,row)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )
            flxz(i,j,k,row) = half*( ( (u0(i,j,k-1,row)**2)/u0(i,j,k-1,ro) ) &
            &               +        ( (u0(i,j,k,row)**2)/u0(i,j,k,ro) )   ) &
            &               + half*( u0(i,j,k-1,pre) + u0(i,j,k,pre) )             
            flxz(i,j,k,roE) = half*( ( (u0(i,j,k-1,roE)*u0(i,j,k-1,row))/u0(i,j,k-1,ro) ) &
            &               +        ( (u0(i,j,k,roE)*u0(i,j,k,rov))/u0(i,j,k,ro) )   )   &
            &               + half*( (u0(i,j,k-1,pre)*u0(i,j,k-1,row)/u0(i,j,k-1,ro)) &
            &               +      ( u0(i,j,k,pre)*u0(i,j,k,row)/u0(i,j,k,ro) ) )
            if(isnan(flxz(i,j,k,ro)) .or. isnan(flxz(i,j,k,rou)) .or. &
            &  isnan(flxz(i,j,k,rov)) .or. isnan(flxz(i,j,k,row)) .or. isnan(flxz(i,j,k,roE))) then
                print*,"i = ", i, " Exiting..NaN found in flxz: ", &
                &            flxz(i,j,k,ro), ", ", flxz(i,j,k,rou), &
                &      ", ", flxz(i,j,k,rov), ", ", flxz(i,j,k,row), ", ", flxz(i,j,k,roE), ", ", flxz(i,j,k,pre) 
                call exit(123)
            endif           
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! zero-order extrapolation for end points of fx
      do n = ro,mach
        do j = fz_lo(2), fz_hi(2)
          do i = fz_lo(1), fz_hi(1)
            flxz(i,j,fz_lo(3),n) = flxz(i,j,fz_lo(3)+1,n)
            flxz(i,j,fz_hi(3),n) = flxz(i,j,fz_hi(3)-1,n)                 
          enddo
        enddo
      enddo

end subroutine compute_con_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute diffusion stage fluxes 
subroutine compute_diff_flux( level, nc, dtdx, dtdy, dtdz, lo, hi,  & 
            &                 u0, u0_lo, u0_hi,               &
            &                 vx, vx_lo, vx_hi,               & 
            &                 vy, vy_lo, vy_hi,               &
            &                 vz, vz_lo, vz_hi,               &
            &                 fldx, fdx_lo, fdx_hi,           &
            &                 fldy, fdy_lo, fdy_hi,           &
            &                 fldz, fdz_lo, fdz_hi            )

  integer, intent(in) :: level, nc
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy, dtdz
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: fdx_lo(3), fdx_hi(3)
  integer, intent(in) :: fdy_lo(3), fdy_hi(3)
  integer, intent(in) :: fdz_lo(3), fdz_hi(3)
  real(amrex_real), intent(inout) :: u0(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),ro:mach)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  real(amrex_real), intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  real(amrex_real), intent(inout) :: fldx(fdx_lo(1):fdx_hi(1),fdx_lo(2):fdx_hi(2),  &
    &                                     fdx_lo(3):fdx_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: fldy(fdy_lo(1):fdy_hi(1),fdy_lo(2):fdy_hi(2),  &
    &                                     fdy_lo(3):fdy_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: fldz(fdz_lo(1):fdz_hi(1),fdz_lo(2):fdz_hi(2),  &
    &                                     fdz_lo(3):fdz_hi(3),ro:pre)


  integer :: i, j, k, n
  real(amrex_real) :: dxdt, dydt, dzdt, epsx, epsy, epsz, nux, nuy, nuz
  ! real(amrex_real), parameter :: one3 = 1.d0/3.d0, one6 = 1.d0/6.d0, half = 0.5_amrex_real

  ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  character(len=128) :: fname, dirchar, levchar

  write(levchar,fmt='(i2.2)') level

  dxdt = 1.d0/dtdx
  dydt = 1.d0/dtdy
  dzdt = 1.d0/dtdz

  ! calculate the diffusion fluxes in x-direction (do not include the last cell)
  !$omp parallel do private(n,k,j,i,epsx,nux)
  do n = ro,roE
    do k = fdx_lo(3), fdx_hi(3)
      do j = fdx_lo(2), fdx_hi(2)
        do i = fdx_lo(1)+1, fdx_hi(1)-1
        ! epsx = dtdx*vx(i,j,k)  
        ! write epsx in terms of velocity at the beginning of the timestep
        ! get vx using the conserved variables
          epsx = dtdx*half*((u0(i-1,j,k,rou)/u0(i-1,j,k,ro)) + (u0(i,j,k,rou)/u0(i,j,k,ro)))
          nux  = half*one6 + one3*(epsx**2)
          fldx(i,j,k,n) = nux*dxdt*(u0(i-1,j,k,n) - u0(i,j,k,n))    
          enddo
        enddo
      enddo
  enddo
  !$omp end parallel do

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
  !$omp parallel do private(n,k,j,i,epsy,nuy)
  do n = ro,roE
    do k = fdy_lo(3), fdy_hi(3)
      do j = fdy_lo(2)+1, fdy_hi(2)-1
        do i = fdy_lo(1), fdy_hi(1)
          epsy = dtdy*half*( (u0(i,j-1,k,rov)/u0(i,j-1,k,ro))  &
          &    + (u0(i,j,k,rov)/u0(i,j,k,ro)) )
          nuy  = one6*half + one3*(epsy**2)
          fldy(i,j,k,n) = nuy*dydt*(u0(i,j-1,k,n) - u0(i,j,k,n))
        enddo 
      enddo
    enddo
  enddo
  !$omp end parallel do

  ! zero-order extrapolation for end points of fldy
  do n = ro, roE
    do k = fdy_lo(3), fdy_hi(3)
      do i = fdy_lo(1), fdy_hi(1)
        fldy(i,fdy_lo(2),k,n) = fldy(i,fdy_lo(2)+1,k,n)
        fldy(i,fdy_hi(2),k,n) = fldy(i,fdy_hi(2)-1,k,n)
      enddo
    enddo
  enddo

  ! calculate the diffusion fluxes in z-direction (do not include the end cells)
  !$omp parallel do private(n,k,j,i,epsz,nuz)
  do n = ro,roE
    do k = fdz_lo(3)+1, fdz_hi(3)-1
      do j = fdz_lo(2), fdz_hi(2)
        do i = fdz_lo(1), fdz_hi(1)
          epsz = dtdz*half*( (u0(i,j,k-1,row)/u0(i,j,k-1,ro))  &
          &    + (u0(i,j,k,row)/u0(i,j,k,ro)) )
          nuz  = one6*half + one3*(epsz**2)
          fldz(i,j,k,n) = nuz*dzdt*(u0(i,j,k-1,n) - u0(i,j,k,n))
        enddo 
      enddo
    enddo
  enddo
  !$omp end parallel do

  ! zero-order extrapolation for end points of fldz
  do n = ro, roE
    do j = fdz_lo(2), fdz_hi(2)
      do i = fdz_lo(1), fdz_hi(1)
        fldz(i,j,fdz_lo(3),n) = fldz(i,j,fdz_lo(3)+1,n)
        fldz(i,j,fdz_hi(3),n) = fldz(i,j,fdz_hi(3)-1,n)
      enddo
    enddo
  enddo

end subroutine compute_diff_flux

!------------------------------------------------------------------------------------------
! Subroutine to compute anti-diffusive fluxes 
subroutine compute_ad_flux( level, time, nc, dtdx, dtdy, dtdz, lo, hi,  & 
      &                 uold, u0_lo, u0_hi,     &
      &                 ucx, ucx_lo, ucx_hi,     &
      &                 ucy, ucy_lo, ucy_hi,     &
      &                 ucz, ucz_lo, ucz_hi,     &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 vz, vz_lo, vz_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi,    &
      &                 fltz, ftz_lo, ftz_hi, diff1     )

  integer, intent(in) :: level, nc
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dtdx, dtdy, dtdz, time, diff1
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: ucx_lo(3), ucx_hi(3)
  integer, intent(in) :: ucy_lo(3), ucy_hi(3)
  integer, intent(in) :: ucz_lo(3), ucz_hi(3)
  integer, intent(in) :: vx_lo(3), vx_hi(3)
  integer, intent(in) :: vy_lo(3), vy_hi(3)
  integer, intent(in) :: vz_lo(3), vz_hi(3)
  integer, intent(in) :: ftx_lo(3), ftx_hi(3)
  integer, intent(in) :: fty_lo(3), fty_hi(3)
  integer, intent(in) :: ftz_lo(3), ftz_hi(3)
  real(amrex_real), intent(inout) :: uold(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
    &                                     ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
    &                                     ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucz (ucz_lo(1):ucz_hi(1),        &
    &                                     ucz_lo(2):ucz_hi(2),ucz_lo(3):ucz_hi(3),ro:pre)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
  real(amrex_real), intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
  real(amrex_real), intent(inout) :: fltx(ftx_lo(1):ftx_hi(1),ftx_lo(2):ftx_hi(2),  &
    &                                     ftx_lo(3):ftx_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: flty(fty_lo(1):fty_hi(1),fty_lo(2):fty_hi(2),  &
    &                                     fty_lo(3):fty_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: fltz(ftz_lo(1):ftz_hi(1),ftz_lo(2):ftz_hi(2),  &
    &                                     ftz_lo(3):ftz_hi(3),ro:roE)

  
  ! integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  ! real(amrex_real), parameter :: half = 0.5_amrex_real, one4 = 0.25_amrex_real, &
  ! &                              one6 = 1.0_amrex_real/6.0_amrex_real
  integer :: i, j, k, n
  integer :: ilo, ihi, jlo, jhi, klo, khi
  real(amrex_real) :: nuxx, nuyy, nuzz, epsx, epsy, epsz
  ! Following variables are used to compute the conserved variables at the corners which is required 
  ! for the prelimiting step
  real(amrex_real) :: ut1(ro:roE), ut2(ro:roE)
  real(amrex_real) :: sgn(ro:roE), du(ro:roE), d2u(ro:roE), fltmp(ro:roE)
  ! real(amrex_real) :: ut1(0:nc-1), ut2(0:nc-1)
  ! real(amrex_real) :: sgn(0:nc-1), du(0:nc-1), d2u(0:nc-1), fltmp(0:nc-1)
  character(len=128) :: dirchar, fname

  ! first calculate x-part of anti-diffusive fluxes (these are at faces normal to x eg. (i+1/2,j,k) )
  ! print*,"ftx_lo= ",ftx_lo, ", ucx_lo = ",ucx_lo
  ! print*,"ftx_hi= ",ftx_hi, ", ucx_hi= ",ucx_hi
  do n = ro,roE
    do k = ftx_lo(3), ftx_hi(3)
      do j = ftx_lo(2), ftx_hi(2)
        do i = ftx_lo(1)+1, ftx_hi(1)-1
          epsx = dtdx*vx(i,j,k)
          nuxx = diff1*(one6 - one6*(epsx**2))
          
          ! ! ! ut1 and ut2 store the values at the corners of the face at which f_ad is calculated
          ! ut1(n) = one4*(ucx(i,j,k,n) + ucx(i,j+1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j+1,k,n))
          ! ut2(n) = one4*(ucx(i,j,k,n) + ucx(i,j-1,k,n) + ucx(i-1,j,k,n) + ucx(i-1,j-1,k,n))

          fltx(i,j,k,n) = nuxx*(ucx(i,j,k,n) - ucx(i-1,j,k,n))  &
          &             - half*one6*(uold(i,j,k,n) - uold(i-1,j,k,n))

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

    do j = ftx_lo(2), ftx_hi(2)
      do i = ftx_lo(1), ftx_hi(1)
        fltx(i,j,ftx_lo(3),n) = fltx(i,j,ftx_lo(3)+1,n)
        fltx(i,j,ftx_hi(3),n) = fltx(i,j,ftx_hi(3)-1,n)        
      enddo
    enddo
  enddo
  
! calculate y-part of anti-diffusive fluxes (these are at faces normal to y eg. (i,j+1/2,k) )
  do n = ro, roE
  ! do n = 0, nc-1
    do k = fty_lo(3),fty_hi(3)
      do j = fty_lo(2)+1,fty_hi(2)-1
        do i = fty_lo(1),fty_hi(1)
          epsy = dtdy*vy(i,j,k)
          nuyy = diff1*(one6 - one6*(epsy**2))
          
          flty(i,j,k,n) = nuyy*(ucy(i,j,k,n) - ucy(i,j-1,k,n))  &
          &             - half*one6*(uold(i,j,k,n) - uold(i,j-1,k,n))

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

    do j = fty_lo(2), fty_hi(2)
      do i = fty_lo(1), fty_hi(1)
        flty(i,j,fty_lo(3),n) = flty(i,j,fty_lo(3)+1,n)
        flty(i,j,fty_hi(3),n) = flty(i,j,fty_hi(3)-1,n)      
      enddo
    enddo 
  enddo

! calculate z-part of anti-diffusive fluxes (these are at faces normal to z eg. (i,j,k+1/2) )
  do n = ro, roE
  ! do n = 0, nc-1
    do k = ftz_lo(3)+1,ftz_hi(3)-1
      do j = ftz_lo(2),ftz_hi(2)
        do i = ftz_lo(1),ftz_hi(1)
          epsz = dtdz*vz(i,j,k)
          nuzz = diff1*(one6 - one6*(epsz**2))
          
          fltz(i,j,k,n) = nuzz*(ucz(i,j,k,n) - ucz(i,j,k-1,n))  &
          &             - half*one6*(uold(i,j,k,n) - uold(i,j,k-1,n))
        enddo
      enddo
    enddo
    ! extrapolate to the end points (zero-order extrapolation)
    do k = ftz_lo(3), ftz_hi(3)
      do j = ftz_lo(2), ftz_hi(2)
        fltz(ftz_lo(1),j,k,n) = fltz(ftz_lo(1)+1,j,k,n)
        fltz(ftz_hi(1),j,k,n) = fltz(ftz_hi(1)-1,j,k,n)
      enddo
      do i = ftz_lo(1), ftz_hi(1)
        fltz(i,ftz_lo(2),k,n) = fltz(i,ftz_lo(2)+1,k,n)
        fltz(i,ftz_hi(2),k,n) = fltz(i,ftz_hi(2)-1,k,n)
      enddo
    enddo

    do j = ftz_lo(2), ftz_hi(2)
      do i = ftz_lo(1), ftz_hi(1)
        fltz(i,j,ftz_lo(3),n) = fltz(i,j,ftz_lo(3)+1,n)
        fltz(i,j,ftz_hi(3),n) = fltz(i,j,ftz_hi(3)-1,n)      
      enddo
    enddo 
  enddo

end subroutine compute_ad_flux

!------------------------------------------------------------------------------------------
! Subroutine to calculate signed quantities and prelimit anti-diffusive fluxes in each direction 

subroutine prelimit_ad_flux( level, nc, lo, hi,  & 
      &                 ucx, ucx_lo, ucx_hi,     &
      &                 ucy, ucy_lo, ucy_hi,     &
      &                 ucz, ucz_lo, ucz_hi,     &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi,    &
      &                 fltz, ftz_lo, ftz_hi     )

  integer, intent(in) :: level, nc
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: ucx_lo(3), ucx_hi(3)
  integer, intent(in) :: ucy_lo(3), ucy_hi(3)
  integer, intent(in) :: ucz_lo(3), ucz_hi(3)
  integer, intent(in) :: ftx_lo(3), ftx_hi(3)
  integer, intent(in) :: fty_lo(3), fty_hi(3)
  integer, intent(in) :: ftz_lo(3), ftz_hi(3)
  real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
    &                                     ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
    &                                     ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: ucz (ucz_lo(1):ucz_hi(1),        &
    &                                     ucz_lo(2):ucz_hi(2),ucz_lo(3):ucz_hi(3),ro:pre)
  real(amrex_real), intent(inout) :: fltx(ftx_lo(1):ftx_hi(1),ftx_lo(2):ftx_hi(2),  &
    &                                     ftx_lo(3):ftx_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: flty(fty_lo(1):fty_hi(1),fty_lo(2):fty_hi(2),  &
    &                                     fty_lo(3):fty_hi(3),ro:roE)
  real(amrex_real), intent(inout) :: fltz(ftz_lo(1):ftz_hi(1),ftz_lo(2):ftz_hi(2),  &
    &                                     ftz_lo(3):ftz_hi(3),ro:roE)

  integer :: i, j, k, n
  integer :: ilo, ihi, jlo, jhi, klo, khi
  real(amrex_real) :: nuxx, nuyy, nuzz, epsx, epsy, epsz
  ! Following variables are used to compute the conserved variables at the corners which is required 
  ! for the prelimiting step
  real(amrex_real) :: ut1(ro:roE), ut2(ro:roE)
  real(amrex_real) :: sgn(ro:roE), du(ro:roE), d2u(ro:roE), fltmp(ro:roE)

! Perform the prelimiting step for flux correction (for x-direction fluxes)
if(level == 0) then
  ilo = lo(1)-1; ihi = hi(1)+1
  jlo = lo(2); jhi = hi(2)
  klo = lo(3); khi = hi(3)
else
  ilo = lo(1)-2; ihi = hi(1)+3
  jlo = lo(2)-3; jhi = hi(2)+3
  klo = lo(3)-3; khi = hi(3)+3
endif
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        fltmp(n) = abs(fltx(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, ucx(i,j,k,n) - ucx(i-1,j,k,n))
        du(n)    = ucx(i-1,j,k,n) - ucx(i-2,j,k,n)
        d2u(n)   = ucx(i+1,j,k,n) - ucx(i,j,k,n) 
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

    do j = ftx_lo(2), ftx_hi(2)
      do i = ftx_lo(1), ftx_hi(1)
        fltx(i,j,ftx_lo(3),n) = fltx(i,j,ftx_lo(3)+1,n)
        fltx(i,j,ftx_hi(3),n) = fltx(i,j,ftx_hi(3)-1,n)        
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
  klo = lo(3)-3; khi = hi(3)+3
endif
! Perform the prelimiting step for flux correction (for y-direction fluxes)
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        fltmp(n) = abs(flty(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, ucy(i,j,k,n) - ucy(i,j-1,k,n))
        du(n)    = ucy(i,j-1,k,n) - ucy(i,j-2,k,n)
        d2u(n)   = ucy(i,j+1,k,n) - ucy(i,j,k,n) 
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
        flty(i,fty_hi(2)-1,k,n) = flty(i,fty_hi(2)-2,k,n)
        flty(i,fty_hi(2),k,n) = flty(i,fty_hi(2)-1,k,n)
      enddo
    enddo
    do j = fty_lo(2), fty_hi(2)
      do i = fty_lo(1), fty_hi(1)
        flty(i,j,fty_lo(3),n) = flty(i,j,fty_lo(3)+1,n)
        flty(i,j,fty_hi(3),n) = flty(i,j,fty_hi(3)-1,n)        
      enddo
    enddo
  endif
enddo
! print*,"prelimited y-antidiffusive fluxes"

if(level == 0) then
  ilo = lo(1); ihi = hi(1)
  jlo = lo(2)-1; jhi = hi(2)+1
  klo = lo(3); khi = hi(3)
else
  ilo = lo(1)-3; ihi = hi(1)+3
  jlo = lo(2)-2; jhi = hi(2)+3
  klo = lo(3)-3; khi = hi(3)+3
endif
! Perform the prelimiting step for flux correction (for y-direction fluxes)
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        fltmp(n) = abs(flty(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, ucy(i,j,k,n) - ucy(i,j-1,k,n))
        du(n)    = ucy(i,j-1,k,n) - ucy(i,j-2,k,n)
        d2u(n)   = ucy(i,j+1,k,n) - ucy(i,j,k,n) 
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
        flty(i,fty_hi(2)-1,k,n) = flty(i,fty_hi(2)-2,k,n)
        flty(i,fty_hi(2),k,n) = flty(i,fty_hi(2)-1,k,n) 
      enddo
    enddo
    do j = fty_lo(2), fty_hi(2)
      do i = fty_lo(1), fty_hi(1)
        flty(i,j,fty_lo(3),n) = flty(i,j,fty_lo(3)+1,n)
        flty(i,j,fty_hi(3),n) = flty(i,j,fty_hi(3)-1,n)        
      enddo
    enddo
  endif
enddo
! print*,"prelimited y-antidiffusive fluxes"

if(level == 0) then
  ilo = lo(1); ihi = hi(1)
  jlo = lo(2); jhi = hi(2)
  klo = lo(3)-1; khi = hi(3)+1
else
  ilo = lo(1)-3; ihi = hi(1)+3
  jlo = lo(2)-3; jhi = hi(2)+3
  klo = lo(3)-2; khi = hi(3)+3
endif
! Perform the prelimiting step for flux correction (for y-direction fluxes)
do n = ro, roE
  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        fltmp(n) = abs(fltz(i,j,k,n))
        sgn(n)   = sign(1.0_amrex_real, ucz(i,j,k,n) - ucz(i,j,k-1,n))
        du(n)    = ucz(i,j,k-1,n) - ucz(i,j,k-2,n)
        d2u(n)   = ucz(i,j,k+1,n) - ucz(i,j,k,n) 
        fltz(i,j,k,n) = sgn(n)*max( 0.0_amrex_real, min(fltmp(n), sgn(n)*du(n), sgn(n)*d2u(n)) )  
      enddo
    enddo
  enddo
! extrapolate to the end points (zero-order extrapolation)
  if(level > 0) then
    do k = ftz_lo(3), ftz_hi(3)
      do j = ftz_lo(2), ftz_hi(2)
        fltz(ftz_lo(1),j,k,n) = fltz(ftz_lo(1)+1,j,k,n)
        fltz(ftz_hi(1),j,k,n) = fltz(ftz_hi(1)-1,j,k,n)
      enddo
      do i = ftz_lo(1), ftz_hi(1)
        fltz(i,ftz_lo(2),k,n) = fltz(i,ftz_lo(2)+1,k,n)
        fltz(i,ftz_hi(2),k,n) = fltz(i,ftz_hi(2)-1,k,n)
      enddo
    enddo
    do j = ftz_lo(2), ftz_hi(2)
      do i = ftz_lo(1), ftz_hi(1)
        fltz(i,j,ftz_lo(3)+1,n) = fltz(i,j,ftz_lo(3)+2,n)
        fltz(i,j,ftz_lo(3),n) = fltz(i,j,ftz_lo(3)+1,n)
        fltz(i,j,ftz_hi(3)-1,n) = fltz(i,j,ftz_hi(3)-2,n)
        fltz(i,j,ftz_hi(3),n) = fltz(i,j,ftz_hi(3)-1,n)        
      enddo
    enddo
  endif
enddo
! print*,"prelimited z-antidiffusive fluxes"

end subroutine prelimit_ad_flux


end module compute_flux_module
