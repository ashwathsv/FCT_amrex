
subroutine advect(level, time, rk, rk_max, fct_step, nc, lo, hi, &
 &            uold, u0_lo, u0_hi, &
 &            ucx, ucx_lo, ucx_hi, &
 &            ucy, ucy_lo, ucy_hi, &
 &            ucz, ucz_lo, ucz_hi, &
 &            uout, uo_lo, uo_hi, &
 &            vx  , vx_lo, vx_hi, &
 &            vy  , vy_lo, vy_hi, &
 &            vz  , vz_lo, vz_hi, &
 &            flxx, fx_lo, fx_hi, &
 &            flxy, fy_lo, fy_hi, &
 &            flxz, fz_lo, fz_hi, &
 &            dx,dt,diff1,pmin,romin) bind(C, name="advect")

use amrex_fort_module, only : amrex_real
use amrex_mempool_module, only : bl_allocate, bl_deallocate
! use compute_flux_module, only : compute_con_flux
use LCPFCT_module, only : LCPFCT3D

implicit none

integer, intent(in) :: level, rk, rk_max, nc, fct_step
integer, intent(in) :: lo(3), hi(3)
real(amrex_real), intent(in) :: dx(3), dt, time, diff1, pmin, romin
integer, intent(in) :: u0_lo(3), u0_hi(3)
integer, intent(in) :: ucx_lo(3), ucx_hi(3)
integer, intent(in) :: ucy_lo(3), ucy_hi(3)
integer, intent(in) :: ucz_lo(3), ucz_hi(3)
integer, intent(in) :: uo_lo(3), uo_hi(3)
integer, intent(in) :: vx_lo(3), vx_hi(3)
integer, intent(in) :: vy_lo(3), vy_hi(3)
integer, intent(in) :: vz_lo(3), vz_hi(3)
integer, intent(in) :: fx_lo(3), fx_hi(3)
integer, intent(in) :: fy_lo(3), fy_hi(3)
integer, intent(in) :: fz_lo(3), fz_hi(3)
real(amrex_real), intent(inout) :: uold(u0_lo(1):u0_hi(1),u0_lo(2):u0_hi(2),u0_lo(3):u0_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: ucx (ucx_lo(1):ucx_hi(1),        &
  &                                       ucx_lo(2):ucx_hi(2),ucx_lo(3):ucx_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: ucy (ucy_lo(1):ucy_hi(1),        &
  &                                       ucy_lo(2):ucy_hi(2),ucy_lo(3):ucy_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: ucz (ucz_lo(1):ucz_hi(1),        &
  &                                       ucz_lo(2):ucz_hi(2),ucz_lo(3):ucz_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: uout(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),0:nc-1)
real(amrex_real), intent(in   ) :: vx  (vx_lo(1):vx_hi(1),vx_lo(2):vx_hi(2),vx_lo(3):vx_hi(3))
real(amrex_real), intent(in   ) :: vy  (vy_lo(1):vy_hi(1),vy_lo(2):vy_hi(2),vy_lo(3):vy_hi(3))
real(amrex_real), intent(in   ) :: vz  (vz_lo(1):vz_hi(1),vz_lo(2):vz_hi(2),vz_lo(3):vz_hi(3))
real(amrex_real), intent(inout) :: flxx(fx_lo(1):fx_hi(1),fx_lo(2):fx_hi(2),fx_lo(3):fx_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: flxy(fy_lo(1):fy_hi(1),fy_lo(2):fy_hi(2),fy_lo(3):fy_hi(3),0:nc-1)
real(amrex_real), intent(inout) :: flxz(fz_lo(1):fz_hi(1),fz_lo(2):fz_hi(2),fz_lo(3):fz_hi(3),0:nc-1)

real(amrex_real) :: umax, vmax, wmax
integer, parameter :: ro = 0, rou = 1, rov = 2, row = 3, roE = 4, pre = 5

! Some compiler may not support 'contiguous'.  Remove it in that case.
! real(amrex_real), dimension(:,:,:,:), pointer, contiguous :: src

! print*,"rk= ",rk,", fct_step= ", fct_step, "nc= ",nc
! print*,"fx_lo= ",fx_lo, ", fx_hi= ",fx_hi
! print*,"fy_lo= ",fy_lo, ", fy_hi= ",fy_hi

umax = maxval(abs(vx))
vmax = maxval(abs(vy))
wmax = maxval(abs(vz))
if ( umax*dt .ge. dx(1) .or. &
     vmax*dt .ge. dx(2) .or. &
     wmax*dt .ge. dx(3) ) then
   print *, "umax = ", umax, ", vmax = ", vmax, ", wmax = ", wmax, ", dt = ", dt, ", dx = ", dx
   call bl_error("CFL violation. Use smaller adv.cfl.")
end if

! first get convective fluxes
    call LCPFCT3D(level, time, fct_step, rk, rk_max,  &
    &           nc, lo, hi,        &
    &           uold, u0_lo, u0_hi,      &
    &           ucx, ucx_lo, ucx_hi,     &
    &           ucy, ucy_lo, ucy_hi,     &
    &           ucz, ucz_lo, ucz_hi,     &
    &           uout, uo_lo, uo_hi,      &
    &           vx,   vx_lo, vx_hi,      &
    &           vy,   vy_lo, vy_hi,      &
    &           vz,   vz_lo, vz_hi,      &
    &           flxx, fx_lo, fx_hi,      &
    &           flxy, fy_lo, fy_hi,      &
    &           flxz, fz_lo, fz_hi,      &
    &           dx, dt, diff1, pmin, romin                   )

end subroutine advect
