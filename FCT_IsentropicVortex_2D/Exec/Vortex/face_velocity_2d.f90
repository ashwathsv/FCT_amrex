!-------------------------------------------------------------------
! Subroutine to calculate face velocities 
subroutine get_face_velocity(level, time, lo, hi, &
     &     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     &     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
     &     dx, prob_lo, state, state_lo, state_hi, &
     &     nc) bind(C, name="get_face_velocity")

  use amrex_fort_module, only : amrex_real
  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: level
  real(amrex_real), intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  integer, intent(in) :: nc
  integer, intent(in) :: state_lo(3), state_hi(3)
  real(amrex_real), intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  real(amrex_real), intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)
  real(amrex_real), intent(in)  :: dx(2), prob_lo(2)
  real(amrex_real), intent(in)  :: state(state_lo(1):state_hi(1), &
  &                                    state_lo(2):state_hi(2), &
  &                                    state_lo(3):state_hi(3), &
  &                                    0:nc-1) 

  integer :: i, j, k

  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  ! calculate vx for all cells (including ghost) except the end face cells
    do k = lo(3), hi(3)
      do j = vx_l2, vx_h2
        do i = vx_l1 + 1, vx_h1 - 1
          vx(i,j) = 0.5_amrex_real*((state(i-1,j,k,rou)/state(i-1,j,k,ro)) + &
          &         (state(i,j,k,rou)/state(i,j,k,ro))) 
        enddo
      enddo
    enddo
    !  the end face values are the same as the values before/after it
    do j = vx_l2, vx_h2
      vx(vx_l1,j) = vx(vx_l1+1,j)
      vx(vx_h1,j) = vx(vx_h1-1,j)
    enddo

  ! calculate vy for all cells (including ghost) except the end face cells
    do k = lo(3), hi(3)
      do j = vy_l2 + 1, vy_h2 - 1
        do i = vy_l1, vy_h1
          vy(i,j) = 0.5_amrex_real*((state(i,j-1,k,rov)/state(i,j-1,k,ro)) + &
          &         (state(i,j,k,rov)/state(i,j,k,ro)))
        enddo
      enddo
    enddo

    !  the end face values are the same as the values before/after it
    do i = vy_l1, vy_h1
      vy(i,vy_l2) = vy(i,vy_l2+1)
      vy(i,vy_h2) = vy(i,vy_h2-1)
    enddo
  
end subroutine get_face_velocity

!-------------------------------------------------------------------
! Subroutine to calculate face velocities to get timestep 

subroutine get_face_velocity_dt(level, time, &
     vx, vx_l1, vx_l2, vx_h1, vx_h2, &
     vy, vy_l1, vy_l2, vy_h1, vy_h2, &
&    dx, prob_lo, umax, vmax, &
&    nc, phi, phi_lo, phi_hi, lo, hi) bind(C, name="get_face_velocity_dt")

  use amrex_fort_module, only : amrex_real
  use amrex_mempool_module, only : bl_allocate, bl_deallocate

  implicit none

  integer, intent(in) :: level, nc
  real(amrex_real), intent(in) :: time
  integer, intent(in) :: vx_l1, vx_l2, vx_h1, vx_h2
  integer, intent(in) :: vy_l1, vy_l2, vy_h1, vy_h2
  integer, intent(in) :: phi_lo(3), phi_hi(3), lo(3), hi(3)
  real(amrex_real), intent(out) :: vx(vx_l1:vx_h1,vx_l2:vx_h2)
  real(amrex_real), intent(out) :: vy(vy_l1:vy_h1,vy_l2:vy_h2)
  real(amrex_real), intent(in)  :: phi(phi_lo(1):phi_hi(1), &
  &                                    phi_lo(2):phi_hi(2), &
  &                                    phi_lo(3):phi_hi(3), &
  &                                    0:nc-1) 
  real(amrex_real), intent(in) :: dx(2), prob_lo(2)
  real(amrex_real), intent(inout) :: umax, vmax

  integer :: i, j, k
  real(amrex_real) :: ss1, ss2, v1, v2
  real(amrex_real), parameter :: M_PI = 3.141592653589793238462643383279502884197d0
  real(amrex_real), parameter :: gamma = 1.4d0
  real(amrex_real), parameter :: c1 = 1.0_amrex_real/(gamma-1)
  real(amrex_real) :: utemp, vtemp

  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

    k = phi_lo(3)
    ! Calculate vx
      do j = vx_l2, vx_h2
        do i = vx_l1, vx_h1
          ! calculate vx at each face
          v1 = phi(i-1,j,k,rou)/phi(i-1,j,k,ro)
          v2 = phi(i,j,k,rou)/phi(i,j,k,ro)
          vx(i,j) = 0.5_amrex_real*(v1 + v2)
            if(isnan(vx(i,j))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in ucx (convection update): ", phi(i,j,k,ro), ", ", phi(i,j,k,rou), &
                &      ", ", phi(i,j,k,rov), ", ", phi(i,j,k,roE), ", ", phi(i,j,k,pre) 
                call exit(123)
            endif           
        enddo
      enddo
      ! Calculate vy
      do j = vy_l2, vy_h2
        do i = vy_l1, vy_h1
          ! calculate vy at each face
          v1 = phi(i,j-1,k,rov)/phi(i,j-1,k,ro)
          v2 = phi(i,j,k,rov)/phi(i,j,k,ro)
          vy(i,j) = 0.5_amrex_real*(v1 + v2)
            if(isnan(vy(i,j))) then
                print*,"location = (", i, ", ",j,"), Exiting..NaN found in ucx (convection update): ", phi(i,j,k,ro), ", ", phi(i,j,k,rou), &
                &      ", ", phi(i,j,k,rov), ", ", phi(i,j,k,roE), ", ", phi(i,j,k,pre) 
                call exit(123)
            endif 
        enddo
      enddo 

    ! calculate maximum velcoity (x and y-components) for stability
      do k = lo(3),hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! get sound speed at the faces of each cell
            ss1 = 0.5_amrex_real*(sqrt(gamma*phi(i,j,k,pre)/phi(i,j,k,ro)) &
              & + sqrt(gamma*phi(i-1,j,k,pre)/phi(i-1,j,k,ro)))
            utemp = (abs(vx(i,j)) + ss1) 
            if (utemp > umax) then
              umax = utemp
            endif          
          enddo
        enddo
      enddo

      do k = lo(3),hi(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            ! get sound speed at the faces of each cell
            ss1 = 0.5_amrex_real*(sqrt(gamma*phi(i,j,k,pre)/phi(i,j,k,ro)) &
              & + sqrt(gamma*phi(i,j-1,k,pre)/phi(i,j-1,k,ro)))
            vtemp = (abs(vy(i,j)) + ss1) 
            if (vtemp > vmax) then
              vmax = vtemp
            endif          
          enddo
        enddo
      enddo   
  
end subroutine get_face_velocity_dt
