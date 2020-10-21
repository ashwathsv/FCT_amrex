module LCPFCT_module

  use amrex_fort_module, only : amrex_real

  implicit none
  
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4, mach = 5
  real(amrex_real), parameter :: one3 = 1.d0/3.d0, one6 = 1.d0/6.d0, gma = 1.4_amrex_real, half = 0.5_amrex_real
  
  private

  public :: LCPFCT2D

  contains

  !----------------------------------------------------------------------------------------
  ! Subroutine to carry out multidimensional FCT
  subroutine LCPFCT2D(level, time, fct_step, rk, rk_max,  &
    &           nc, lo, hi,        &
    &           uold, u0_lo, u0_hi,      &
    &           ucx, ucx_lo, ucx_hi,     &
    &           ucy, ucy_lo, ucy_hi,     &
    &           uout, uo_lo, uo_hi,      &
    &           vx,   vx_lo, vx_hi,      &
    &           vy,   vy_lo, vy_hi,      &
    &           flxx, fx_lo, fx_hi,      &
    &           flxy, fy_lo, fy_hi,      &
    &           dx, dt                   ) 

  use amrex_fort_module, only : amrex_real
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use compute_flux_module, only : compute_con_flux, compute_diff_flux, compute_ad_flux, compute_source_flux

  integer, intent(in) :: level, nc, fct_step, rk, rk_max
  integer, intent(in) :: lo(3), hi(3)
  real(amrex_real), intent(in) :: dx(2), dt, time
  integer, intent(in) :: u0_lo(3), u0_hi(3)
  integer, intent(in) :: ucx_lo(3), ucx_hi(3)
  integer, intent(in) :: ucy_lo(3), ucy_hi(3)
  integer, intent(in) :: uo_lo(3), uo_hi(3)
  integer, intent(in) :: vx_lo(2), vx_hi(2)
  integer, intent(in) :: vy_lo(2), vy_hi(2)
  integer, intent(in) :: fx_lo(3), fx_hi(3)
  integer, intent(in) :: fy_lo(3), fy_hi(3)
  real(amrex_real), intent(inout) :: uold(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
    &                                       ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
    &                                       ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),0:nc-1)
  real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
  real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
  real(amrex_real), intent(inout) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:nc-1)
  real(amrex_real), intent(inout) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:nc-1)


  integer :: i, j, k, n
  integer :: ilo, ihi, jlo, jhi, klo, khi
  real(amrex_real) :: dtdx, dtdy, coeff, dxdt, dydt, velmod, ss
  integer :: ftx_lo(3), ftx_hi(3)
  integer :: fty_lo(3), fty_hi(3)
  integer :: ut_lo(3), ut_hi(3) 

  real(amrex_real), dimension(:,:,:,:), pointer, contiguous :: fltx, flty, utemp, frin, frout
  real(amrex_real), dimension(:), pointer, contiguous :: umin, umax, flin, flout, temp

  character(len=128) :: fname, levchar, rkchar
  write(levchar,fmt='(i2.2)') level
  write(rkchar,fmt='(i2.2)') rk
  if(rk_max == 4) then
    if(rk == 1) then
      coeff = 0.25_amrex_real
    else if(rk==2) then
      coeff = one3
    else if (rk == 3) then
      coeff = 0.5_amrex_real
    else
      coeff = 1.0_amrex_real
    endif
  else if(rk_max == 2) then
    if(rk == 1) then
      coeff = 0.5_amrex_real
    else
      coeff = 1.0_amrex_real
    endif
  else 
    coeff = 1.0_amrex_real
  endif

  dtdx = coeff*dt/dx(1)
  dtdy = coeff*dt/dx(2)
  dxdt = 1.0_amrex_real/dtdx
  dydt = 1.0_amrex_real/dtdy

  if(fct_step == 1) then
    call compute_con_flux(  level, nc, lo, hi,  & 
        &                 uold, u0_lo, u0_hi,       &
        &                 flxx, fx_lo, fx_hi,       &
        &                 flxy, fy_lo, fy_hi  )
    ! print*,"rk= ",rk,",dx(1)= ",dx(1),", from LCPFCT2D, fx_lo= ",fx_lo,", fx_hi= ",fx_hi
    ! Predictor step of FCT
    ! compute convected x and y values of the conserved quantities
    if(level == 0) then
      ilo = lo(1);  ihi = hi(1)
      jlo = lo(2);  jhi = hi(2)
      klo = lo(3);  khi = hi(3)
    else
      ilo = lo(1) - 3;  ihi = hi(1) + 3
      jlo = lo(2) - 3;  jhi = hi(2) + 3
      klo = lo(3);  khi = hi(3)
    endif

    ! deciding which quantity is to be used for calculating fluxes in different rk steps
    if(rk == 1) then
      ut_lo = u0_lo;  ut_hi = u0_hi
      call bl_allocate(utemp,ut_lo(1),ut_hi(1),ut_lo(2),ut_hi(2),ut_lo(3),ut_hi(3),ro,pre)
      utemp = uold(:,:,:,ro:pre)
    else
      ut_lo = uo_lo;  ut_hi = uo_hi
      call bl_allocate(utemp,ut_lo(1),ut_hi(1),ut_lo(2),ut_hi(2),ut_lo(3),ut_hi(3),ro,pre)
      utemp = uout(:,:,:,ro:pre)
    endif

    do n = ro,roE ! do not update pressure here
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            ucx(i,j,k,n)  = uold(i,j,k,n) - dtdx*(flxx(i+1,j,k,n) - flxx(i,j,k,n))
            ucy(i,j,k,n)  = uold(i,j,k,n) - dtdy*(flxy(i,j+1,k,n) - flxy(i,j,k,n))
            uout(i,j,k,n) = uold(i,j,k,n) - dtdx*(flxx(i+1,j,k,n) - flxx(i,j,k,n)) &
            &             - dtdy*(flxy(i,j+1,k,n) - flxy(i,j,k,n))
            if(isnan(ucx(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in ucx (convection update): ", ucx(i,j,k,ro), ", ", ucx(i,j,k,rou), &
                &      ", ", ucx(i,j,k,rov), ", ", ucx(i,j,k,roE), ", ", ucx(i,j,k,pre) 
                call exit(123)
            endif 
            if(isnan(ucy(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in ucy (convection update): ", ucy(i,j,k,ro), ", ", ucy(i,j,k,rou), &
                &      ", ", ucy(i,j,k,rov), ", ", ucy(i,j,k,roE), ", ", ucy(i,j,k,pre) 
                call exit(123)
            endif 
            if(isnan(uout(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in uout (convection update): ", uout(i,j,k,ro), ", ", uout(i,j,k,rou), &
                &      ", ", uout(i,j,k,rov), ", ", uout(i,j,k,roE), ", ", uout(i,j,k,pre) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo
    enddo

    if(level > 0) then
      ! zero order extrapolation for end points
      ! do n = 0, nc-3
      do n = ro, pre
        do k = klo, khi
          do j = jlo, jhi
            ucx(ilo-1,j,k,n)  = ucx(ilo,j,k,n);  ucx(ihi+1,j,k,n)  = ucx(ihi,j,k,n)
            ucy(ilo-1,j,k,n)  = ucy(ilo,j,k,n);  ucy(ihi+1,j,k,n)  = ucy(ihi,j,k,n)
            uout(ilo-1,j,k,n) = uout(ilo,j,k,n); uout(ihi+1,j,k,n) = uout(ihi,j,k,n)
          enddo
          do i = ilo-1,ihi+1           
            ucx(i,jlo-1,k,n)  = ucx(i,jlo,k,n);  ucx(i,jhi+1,k,n)  = ucx(i,jhi,k,n)
            ucy(i,jlo-1,k,n)  = ucy(i,jlo,k,n);  ucy(i,jhi+1,k,n)  = ucy(i,jhi,k,n)
            uout(i,jlo-1,k,n) = uout(i,jlo,k,n); uout(i,jhi+1,k,n) = uout(i,jhi,k,n)
          enddo
        enddo
      enddo
    endif

    ! convective step good for 2 level runs (same results in x, y propagation)
    ! allocate arrays for diffusion stage
    ftx_lo = fx_lo; ftx_hi = fx_hi
    fty_lo = fy_lo; fty_hi = fy_hi

    call bl_allocate(fltx,ftx_lo(1),ftx_hi(1),ftx_lo(2),ftx_hi(2),ftx_lo(3),ftx_hi(3),ro,roE)
    call bl_allocate(flty,fty_lo(1),fty_hi(1),fty_lo(2),fty_hi(2),fty_lo(3),fty_hi(3),ro,roE)

    call compute_diff_flux(  level, nc, dtdx, dtdy, lo, hi,  & 
      &                 uold, u0_lo, u0_hi,      &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi     )

    if(rk == rk_max) then
      do n = ro,roE
        flxx(:,:,:,n)  = flxx(:,:,:,n) + fltx(:,:,:,n)
        flxy(:,:,:,n)  = flxy(:,:,:,n) + flty(:,:,:,n) 
      enddo
    endif

    do n = ro,roE  
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            uout(i,j,k,n) = uout(i,j,k,n) - dtdx*(fltx(i+1,j,k,n) - fltx(i,j,k,n)) &
            &             - dtdy*(flty(i,j+1,k,n) - flty(i,j,k,n))
            if(isnan(uout(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in uout (diffusion update): ", uout(i,j,k,ro), ", ", uout(i,j,k,rou), &
                &      ", ", uout(i,j,k,rov), ", ", uout(i,j,k,roE), ", ", uout(i,j,k,pre) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo
    enddo

    if(level > 0) then
      ! zero order extrapolation for end points
      do n = ro,roE
      ! do n = 0,nc-1
        do k = klo, khi
          do j = jlo, jhi
            uout(ilo-1,j,k,n) = uout(ilo,j,k,n); uout(ihi+1,j,k,n) = uout(ihi,j,k,n);
          enddo
          do i = ilo-1,ihi+1           
            uout(i,jlo-1,k,n) = uout(i,jlo,k,n); uout(i,jhi+1,k,n) = uout(i,jhi,k,n);
          enddo
        enddo
      enddo
    endif

    ! diffusion step also works fine (same solution for x, y propagation)

    ! compute source terms (these variables are stored again in fltx and flty)
    call compute_source_flux( level, nc, dtdx, dtdy, lo, hi,  & 
      &                 utemp, ut_lo, ut_hi,      &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi     )

    if(rk == rk_max) then
      do n = ro,roE
        flxx(:,:,:,n)  = flxx(:,:,:,n) + fltx(:,:,:,n)
        flxy(:,:,:,n)  = flxy(:,:,:,n) + flty(:,:,:,n) 
      enddo
    endif

    ! update source terms and store results (ro^l) in uout
    do n = ro,roE
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            uout(i,j,k,n) = uout(i,j,k,n) - dtdx*(fltx(i+1,j,k,n) - fltx(i,j,k,n))  &
            &             - dtdy*(flty(i,j+1,k,n) - flty(i,j,k,n))
            if(isnan(uout(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in uout (source update): ", uout(i,j,k,ro), ", ", uout(i,j,k,rou), &
                &      ", ", uout(i,j,k,rov), ", ", uout(i,j,k,roE), ", ", uout(i,j,k,pre) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo
    enddo

    if(level > 0) then
      ! zero order extrapolation for end points (only for levels other than level 0)
      do n = ro,roE
        do k = klo, khi
          do j = jlo, jhi
            uout(ilo-1,j,k,n) = uout(ilo,j,k,n); uout(ihi+1,j,k,n) = uout(ihi,j,k,n);
          enddo
          do i = ilo-1,ihi+1           
            uout(i,jlo-1,k,n) = uout(i,jlo,k,n); uout(i,jhi+1,k,n) = uout(i,jhi,k,n);
          enddo
        enddo
      enddo
    endif  

    call bl_deallocate(utemp)
    call bl_deallocate(fltx)
    call bl_deallocate(flty)

  else

    ut_lo = uo_lo;  ut_hi = uo_hi
    ftx_lo = fx_lo; ftx_hi = fx_hi
    fty_lo = fy_lo; fty_hi = fy_hi
    ! allocate arrays for anti-diffusion stage
    call bl_allocate(fltx,ftx_lo(1),ftx_hi(1),ftx_lo(2),ftx_hi(2),ftx_lo(3),ftx_hi(3),ro,roE)
    call bl_allocate(flty,fty_lo(1),fty_hi(1),fty_lo(2),fty_hi(2),fty_lo(3),fty_hi(3),ro,roE)
    call bl_allocate(utemp,ut_lo(1),ut_hi(1),ut_lo(2),ut_hi(2),ut_lo(3),ut_hi(3),ro,pre)
    utemp = uout(:,:,:,ro:pre)

    ! compute antidiffusive fluxes (these variables are stored again in fltx and flty) and do the 
    ! prelimiting step
    call compute_ad_flux( level, time, nc, dtdx, dtdy, lo, hi,  & 
      &                 utemp, ut_lo, ut_hi,     &
      &                 ucx, ucx_lo, ucx_hi,     &
      &                 ucy, ucy_lo, ucy_hi,     &
      &                 vx, vx_lo, vx_hi,        & 
      &                 vy, vy_lo, vy_hi,        &
      &                 fltx, ftx_lo, ftx_hi,    &
      &                 flty, fty_lo, fty_hi     )

    if(level == 0) then
      ilo = lo(1)-1; ihi = hi(1)+1
      jlo = lo(2)-1; jhi = hi(2)+1
      klo = lo(3);   khi = hi(3)
    else
      ilo = lo(1)-3; ihi = hi(1)+3
      jlo = lo(2)-3; jhi = hi(2)+3
      klo = lo(3);   khi = hi(3)
    endif

    call bl_allocate(umin,ro,roE)
    call bl_allocate(umax,ro,roE)
    call bl_allocate(frin,ilo,ihi,jlo,jhi,ut_lo(3),ut_hi(3),ro,roE)
    call bl_allocate(frout,ilo,ihi,jlo,jhi,ut_lo(3),ut_hi(3),ro,roE)
    call bl_allocate(flin,ro,roE)
    call bl_allocate(flout,ro,roE)

    ! Flux correction procedure (steps A, C-F in Devore)
    do n = ro,roE
    ! do n = 0,nc-1
      do k = lo(3), hi(3)
        do j = jlo, jhi
          do i = ilo, ihi
            ! Limits for conserved variables
            umin(n) = min(utemp(i,j-1,k,n), utemp(i-1,j,k,n), utemp(i,j,k,n), utemp(i+1,j,k,n), utemp(i,j+1,k,n))
            umax(n) = max(utemp(i,j-1,k,n), utemp(i-1,j,k,n), utemp(i,j,k,n), utemp(i+1,j,k,n), utemp(i,j+1,k,n))

            ! calculate total incoming and outgoing antidiffusive fluxes in each cell
            flin(n) = max(fltx(i,j,k,n),0.d0) - min(fltx(i+1,j,k,n),0.d0) &
            &       + max(flty(i,j,k,n),0.d0) - min(flty(i,j+1,k,n),0.d0)

            flout(n) = max(fltx(i+1,j,k,n),0.d0) - min(fltx(i,j,k,n),0.d0) &
            &        + max(flty(i,j+1,k,n),0.d0) - min(flty(i,j,k,n),0.d0)

            ! calculate fractions of incoming and outgoing fluxes applied to each cell
              frin(i,j,k,n) = (umax(n) - utemp(i,j,k,n))/(1E-16_amrex_real + flin(n))
              frout(i,j,k,n) = (utemp(i,j,k,n) - umin(n))/(1E-16_amrex_real + flout(n))

            if(isnan(frin(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in frin (step C): ", &
                &       "n= ", n, ", flin= ", flin(n), ", frin= ", frin(i,j,k,ro), ", ", frin(i,j,k,rou), &
                &      ", ", frin(i,j,k,rov), ", ", frin(i,j,k,roE) 
                call exit(123)
            endif 
            if(isnan(frout(i,j,k,n))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in frout (step C): ", &
                &       "n= ", n, ", flin= ", flout(n), ", frin= ", frout(i,j,k,ro), ", ", frout(i,j,k,rou), &
                &      ", ", frout(i,j,k,rov), ", ", frout(i,j,k,roE) 
                call exit(123)
            endif 
          enddo
        enddo
      enddo
    enddo

    call bl_allocate(temp,ro,roE)
    ! call bl_allocate(temp,0,nc-1)
    if(level == 0) then
      ilo = lo(1); ihi = hi(1)+1
      jlo = lo(2); jhi = hi(2)
      klo = lo(3);   khi = hi(3)
    else
      ilo = lo(1)-2; ihi = hi(1)+3
      jlo = lo(2)-3; jhi = hi(2)+3
      klo = lo(3);   khi = hi(3)
    endif
    ! calculate the corrected fluxes before updating the conserved variables
    do n = ro,roE
      ! update fluxes at faces whose normals are in x-direction (fltx)
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            temp(n) = fltx(i,j,k,n)
            if(temp(n) >= 0.d0) then
              fltx(i,j,k,n) = temp(n)*min(frout(i-1,j,k,n),frin(i,j,k,n),1.d0)
            else
              fltx(i,j,k,n) = temp(n)*min(frin(i-1,j,k,n),frout(i,j,k,n),1.d0)
            endif
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

    if(level == 0) then
      ilo = lo(1); ihi = hi(1)
      jlo = lo(2); jhi = hi(2)+1
      klo = lo(3);   khi = hi(3)
    else
      ilo = lo(1)-3; ihi = hi(1)+3
      jlo = lo(2)-2; jhi = hi(2)+3
      klo = lo(3);   khi = hi(3)
    endif
      ! update fluxes at faces whose normals are in y-direction (flty)
    do n = ro,roE
      do k = klo, khi
        do j = jlo, jhi
          do i = ilo, ihi
            temp(n) = flty(i,j,k,n)
            if(temp(n) >= 0.d0) then
              flty(i,j,k,n) = temp(n)*min(frout(i,j-1,k,n),frin(i,j,k,n),1.d0)
            else
              flty(i,j,k,n) = temp(n)*min(frin(i,j-1,k,n),frout(i,j,k,n),1.d0)
            endif
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
    endif
  enddo

    ! update conserved variables
    if(level == 0) then
      ilo = lo(1); ihi = hi(1)
      jlo = lo(2); jhi = hi(2)
      klo = lo(3); khi = hi(3)
    else
      ilo = lo(1)-3; ihi = hi(1)+3
      jlo = lo(2)-3; jhi = hi(2)+3
      klo = lo(3);   khi = hi(3)
    endif

    do n = ro,roE
      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            uout(i,j,k,n) = utemp(i,j,k,n) - (fltx(i+1,j,k,n) - fltx(i,j,k,n)) &
            &             - (flty(i,j+1,k,n) - flty(i,j,k,n))
          enddo
        enddo
      enddo
    enddo

    ! update pressure and entropy
    do k = klo,khi
      do j = jlo,jhi
        do i = ilo,ihi
          uout(i,j,k,pre) = (gma-1)*( uout(i,j,k,roE)                       &
          &               -  half*( (uout(i,j,k,rou)**2 + uout(i,j,k,rov)**2)/uout(i,j,k,ro) ) )

          ss = sqrt(gma*uout(i,j,k,pre)/uout(i,j,k,ro))
          velmod = sqrt( (uout(i,j,k,rou)/uout(i,j,k,ro))**2 + (uout(i,j,k,rov)/uout(i,j,k,ro))**2 )
          uout(i,j,k,mach) = velmod/ss
        enddo
      enddo
    enddo
    ! extrapolate to end points (zero-order extrapolation)
    if(level > 0) then
      do n = ro,mach
        do k = klo,khi
          do j = uo_lo(2)+1,uo_hi(2)-1
            uout(uo_lo(1),j,k,n) = uout(uo_lo(1)+1,j,k,n)
            uout(uo_hi(1),j,k,n) = uout(uo_hi(1)-1,j,k,n)
          enddo
          do i = uo_lo(1),uo_hi(1)
            uout(i,uo_lo(2),k,n) = uout(i,uo_lo(2)+1,k,n)
            uout(i,uo_hi(2),k,n) = uout(i,uo_hi(2)-1,k,n)
          enddo
        enddo
      enddo
    endif         

    ! if(ddir == 1) then
      ! print*,"before scaling, end of rk= ",rk,", lev= ,",level,", max(abs(flxy))= ",maxval(abs(flxy(:,:,:,rov)))
      ! print*,"level= ",level,"dx= ",dx
    ! else
      ! print*,"end of rk= ",rk,", lev= ,",level,", max(abs(flxx))= ",maxval(abs(flxx(:,:,:,rou)))
    ! endif

    ! scale fluxes by time and area
    if(rk == rk_max) then
      flxx(:,:,:,pre:mach) = 0.0_amrex_real
      flxy(:,:,:,pre:mach) = 0.0_amrex_real

      do n = ro,roE
        ! scale x-fluxes
        do k = fx_lo(3), fx_hi(3)
          do j = fx_lo(2), fx_hi(2)
            do i = fx_lo(1), fx_hi(1)
              flxx(i,j,k,n) = (flxx(i,j,k,n) + dxdt*fltx(i,j,k,n))*dx(2)*dt
            enddo
          enddo
        enddo
        
        ! scale y-fluxes
        do k = fy_lo(3), fy_hi(3)
          do j = fy_lo(2), fy_hi(2)
            do i = fy_lo(1), fy_hi(1)
              flxy(i,j,k,n) = (flxy(i,j,k,n) + dydt*flty(i,j,k,n))*dx(1)*dt
            enddo
          enddo
        enddo

      enddo
    endif
    
    call bl_deallocate(fltx)
    call bl_deallocate(flty)
    call bl_deallocate(utemp)
    call bl_deallocate(umin)
    call bl_deallocate(umax)
    call bl_deallocate(umin)
    call bl_deallocate(frin)
    call bl_deallocate(frout)
    call bl_deallocate(flin)
    call bl_deallocate(flout)
    call bl_deallocate(temp)

  endif

end subroutine LCPFCT2D

end module LCPFCT_module
