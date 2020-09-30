
subroutine advect(level, time, rk, rk_max, fct_step, ddir, nc, lo, hi, &
 &            uold, u0_lo, u0_hi, &
 &            ucx, ucx_lo, ucx_hi, &
 &            ucy, ucy_lo, ucy_hi, &
 &            uout, uo_lo, uo_hi, &
 &            vx  , vx_lo, vx_hi, &
 &            vy  , vy_lo, vy_hi, &
 &            flxx, fx_lo, fx_hi, &
 &            flxy, fy_lo, fy_hi, &
 &            dx,dt) bind(C, name="advect")

use amrex_fort_module, only : amrex_real
use amrex_mempool_module, only : bl_allocate, bl_deallocate
use compute_flux_module, only : compute_con_flux
use LCPFCT_module, only : LCPFCT2D

implicit none

integer, intent(in) :: level, rk, rk_max, nc, ddir, fct_step
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

real(amrex_real) :: umax, vmax
integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

! Some compiler may not support 'contiguous'.  Remove it in that case.
! real(amrex_real), dimension(:,:,:,:), pointer, contiguous :: src

! print*,"rk= ",rk,", fct_step= ", fct_step, "nc= ",nc
! print*,"fx_lo= ",fx_lo, ", fx_hi= ",fx_hi
! print*,"fy_lo= ",fy_lo, ", fy_hi= ",fy_hi

umax = maxval(abs(vx))
vmax = maxval(abs(vy))
if ( umax*dt .ge. dx(1) .or. &
     vmax*dt .ge. dx(2) ) then
   print *, "umax = ", umax, ", vmax = ", vmax, ", dt = ", dt, ", dx = ", dx
   call bl_error("CFL violation. Use smaller adv.cfl.")
end if

! first get convective fluxes
    call LCPFCT2D(level, time, fct_step, rk, rk_max,  &
    &           ddir, nc, lo, hi,        &
    &           uold, u0_lo, u0_hi,      &
    &           ucx, ucx_lo, ucx_hi,     &
    &           ucy, ucy_lo, ucy_hi,     &
    &           uout, uo_lo, uo_hi,      &
    &           vx,   vx_lo, vx_hi,      &
    &           vy,   vy_lo, vy_hi,      &
    &           flxx, fx_lo, fx_hi,      &
    &           flxy, fy_lo, fy_hi,      &
    &           dx, dt                   )
! ! edge states
! call bl_allocate(phix_1d, glo(1), ghi(1), glo(2), ghi(2))
! call bl_allocate(phiy_1d, glo(1), ghi(1), glo(2), ghi(2))
! call bl_allocate(phix   , glo(1), ghi(1), glo(2), ghi(2))
! call bl_allocate(phiy   , glo(1), ghi(1), glo(2), ghi(2))
! ! slope                                                 
! call bl_allocate(slope  , glo(1), ghi(1), glo(2), ghi(2))

! ! We like to allocate these **pointers** here and then pass them to a function
! ! to remove their pointerness for performance, because normally pointers could
! ! be aliasing.  We need to use pointers instead of allocatable arrays because
! ! we like to use AMReX's bl_allocate to allocate memeory instead of the intrinsic
! ! allocate.  Bl_allocate is much faster than allocate inside OMP.  
! ! Note that one MUST CALL BL_DEALLOCATE.

! ! check if CFL condition is violated.


! ! call a function to compute flux
! call compute_flux_2d(lo, hi, dt, dx, &
!                      uin, ui_lo, ui_hi, &
!                      vx, vx_lo, vx_hi, &
!                      vy, vy_lo, vy_hi, &
!                      flxx, fx_lo, fx_hi, &
!                      flxy, fy_lo, fy_hi, &
!                      phix_1d, phiy_1d, phix, phiy, slope, glo, ghi)

! ! Do a conservative update
! do    j = lo(2),hi(2)
!    do i = lo(1),hi(1)
!       uout(i,j) = uin(i,j) + &
!            ( (flxx(i,j) - flxx(i+1,j)) * dtdx(1) &
!            + (flxy(i,j) - flxy(i,j+1)) * dtdx(2) )
!    enddo
! enddo

! ! Scale by face area in order to correctly reflx
! do    j = lo(2), hi(2)
!    do i = lo(1), hi(1)+1
!       flxx(i,j) = flxx(i,j) * ( dt * dx(2))
!    enddo
! enddo

! ! Scale by face area in order to correctly reflx
! do    j = lo(2), hi(2)+1 
!    do i = lo(1), hi(1)
!       flxy(i,j) = flxy(i,j) * (dt * dx(1))
!    enddo
! enddo

! call bl_deallocate(phix_1d)
! call bl_deallocate(phiy_1d)
! call bl_deallocate(phix)
! call bl_deallocate(phiy)
! call bl_deallocate(slope)

end subroutine advect
!----------------------------------------------------------------------------------------
! Subroutine to carry out multidimensional FCT
! subroutine LCPFCT2D(level, fct_step, coeff,  &
!         &           ddir, nc, lo, hi,        &
!         &           uold, u0_lo, u0_hi,      &
!         &           ucx, ucx_lo, ucx_hi,     &
!         &           ucy, ucy_lo, ucy_hi,     &
!         &           uold, uo_lo, uo_hi,      &
!         &           vx,   vx_lo, vx_hi,      &
!         &           vy,   vy_lo, vy_hi,      &
!         &           flxx, fx_lo, fx_hi,      &
!         &           flxy, fy_lo, fy_hi,      &
!         &           dx, dt                   ) 

!   use amrex_fort_module, only : amrex_real
!   use amrex_mempool_module, only : bl_allocate, bl_deallocate
!   use compute_flux_module, only : compute_flux_2d, compute_con_flux

!   implicit none

!   integer, intent(in) :: level, nc, ddir, fct_step
!   integer, intent(in) :: lo(3), hi(3)
!   real(amrex_real), intent(in) :: dx(2), dt, coeff
!   integer, intent(in) :: u0_lo(3), u0_hi(3)
!   integer, intent(in) :: ucx_lo(3), ucx_hi(3)
!   integer, intent(in) :: ucy_lo(3), ucy_hi(3)
!   integer, intent(in) :: uo_lo(3), uo_hi(3)
!   integer, intent(in) :: vx_lo(2), vx_hi(2)
!   integer, intent(in) :: vy_lo(2), vy_hi(2)
!   integer, intent(in) :: fx_lo(3), fx_hi(3)
!   integer, intent(in) :: fy_lo(3), fy_hi(3)
!   real(amrex_real), intent(inout) :: uold(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
!   real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
!   &                                       ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),0:nc-1)
!   real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
!   &                                       ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),0:nc-1)
!   real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),0:nc-1)
!   real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2))
!   real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2))
!   real(amrex_real), intent(inout) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:nc-1)
!   real(amrex_real), intent(inout) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:nc-1)


!   integer :: i, j, k
!   integer :: ilo, ihi, jlo, jhi, jlo, khi
!   real(amrex_real) :: dtdx, dtdy

!   ! dtdx = coeff*dt/dx(1)
!   ! dtdy = coeff*dt/dx(2)

!   ! if(fct_step == 1) then
!   ! ! Predictor step of FCT

!   ! ! compute convected x and y values of the conserved quantities
!   ! if(level == 0) then
!   !   ilo = lo(1);  ihi = hi(1)
!   !   jlo = lo(2);  jhi = hi(2)
!   !   klo = lo(3);  khi = hi(3)
!   ! else
!   !   ilo = lo(1) - 3;  ihi = hi(1) + 3
!   !   jlo = lo(2) - 3;  jhi = hi(2) + 3
!   !   klo = lo(3);  khi = hi(3)
!   ! endif

! end subroutine LCPFCT2D