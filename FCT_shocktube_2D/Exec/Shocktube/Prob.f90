
subroutine initdata(level, time, nc, domdir, lo, hi, &
     phi, phi_lo, phi_hi, &
     dx, prob_lo, prob_hi, ro2ro1, v2v1, p2p1) bind(C, name="initdata")

  use amrex_fort_module, only : amrex_spacedim, amrex_real

  implicit none
  integer, intent(in) :: level, nc, domdir, lo(3), hi(3), phi_lo(3), phi_hi(3)
  real(amrex_real), intent(in) :: time, ro2ro1, v2v1, p2p1
  real(amrex_real), intent(inout) :: phi(phi_lo(1):phi_hi(1), &
       &                                 phi_lo(2):phi_hi(2), &
       &                                 phi_lo(3):phi_hi(3), &
       &                                 0:nc-1)
  real(amrex_real), intent(in) :: dx(3), prob_lo(3), prob_hi(3)

  integer          :: i,j,k
  real(amrex_real) :: x,y,z,r2,mid
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4
  real(amrex_real), parameter :: gamma = 1.4
  real(amrex_real), parameter :: c1 = 1.0_amrex_real/(gamma-1)
  
  !$omp parallel do private(i,j,k,x,y,z,r2) collapse(2)
  do k=lo(3),hi(3)
     do j=lo(2),hi(2)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1),hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

            if(domdir == 1) then
              mid = 0.5_amrex_real*(prob_lo(1) + prob_hi(1))
              phi(i,j,k,rov) = 0.0_amrex_real
              if(x < mid) then
                phi(i,j,k,ro)  = ro2ro1
                phi(i,j,k,rou) = v2v1 
                phi(i,j,k,pre) = p2p1
              else
                phi(i,j,k,ro)  = 1.0_amrex_real
                phi(i,j,k,rou) = 0.0_amrex_real
                phi(i,j,k,pre) = 1.0_amrex_real
              endif
                phi(i,j,k,roE) = phi(i,j,k,pre)*c1 + &
                &                0.5_amrex_real*((phi(i,j,k,rou)**2)/phi(i,j,k,ro))                
            else if(domdir == 2) then
              mid = 0.5_amrex_real*(prob_lo(2) + prob_hi(2))
              phi(i,j,k,rou) = 0.0_amrex_real
              if(y < mid) then
                phi(i,j,k,ro)  = ro2ro1
                phi(i,j,k,rov) = v2v1
                phi(i,j,k,pre) = p2p1
              else
                phi(i,j,k,ro)  = 1.0_amrex_real
                phi(i,j,k,rov) = 0.0_amrex_real
                phi(i,j,k,pre) = 1.0_amrex_real                
              endif
                phi(i,j,k,roE) = phi(i,j,k,pre)*c1 + &
                &                0.5_amrex_real*((phi(i,j,k,rou)**2)/phi(i,j,k,ro))  
            endif

        end do
     end do
  end do
  !$omp end parallel do

end subroutine initdata
