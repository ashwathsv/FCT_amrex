
subroutine initdata(level, time, nc, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo, prob_hi, ro2, ro1, u, v, p, xw, yw) bind(C, name="initdata")

  use amrex_fort_module, only : amrex_spacedim, amrex_real

  implicit none
  integer, intent(in) :: level, nc, lo(3), hi(3), phi_lo(3), phi_hi(3)
  real(amrex_real), intent(in) :: time, ro2, ro1, u, v, p, xw, yw
  real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3), &
       &                                 0:nc-1)
  real(amrex_real), intent(in) :: dx(3), prob_lo(3), prob_hi(3)

  integer          :: i,j,k
  real(amrex_real) :: x,y,z,mid(2)
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  real(amrex_real), parameter :: gamma = 1.4_amrex_real, half = 0.5_amrex_real
  real(amrex_real), parameter :: c1 = 1.0_amrex_real/(gamma-1)

  mid(1) = 0.5_amrex_real*(prob_lo(1) + prob_hi(1))
  mid(2) = 0.5_amrex_real*(prob_lo(2) + prob_hi(2))
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)
            if(x <= mid(1)+half*xw .and. x >= mid(1)-half*xw .and.  &
            &  y <= mid(2)+half*yw .and. y >= mid(2)-half*yw) then
              phi(i,j,k,ro) = ro2
              phi(i,j,k,rou) = ro2*u
              phi(i,j,k,rov) = ro2*v
            else
              phi(i,j,k,ro) = ro1
              phi(i,j,k,rou) = ro1*u
              phi(i,j,k,rov) = ro1*v
            endif
            phi(i,j,k,pre) = p
            phi(i,j,k,roE) = p*c1 + half*((phi(i,j,k,rou)**2 + phi(i,j,k,rov)**2)/phi(i,j,k,ro))
        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
