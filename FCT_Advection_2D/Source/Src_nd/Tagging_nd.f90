
! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the state
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag        <=  integer tag array
! ::: tag_lo,hi   => index extent of tag array
! ::: state       => state array
! ::: state_lo,hi => index extent of state array
! ::: set         => integer value to tag cell for refinement
! ::: clear       => integer value to untag cell
! ::: lo,hi       => work region we are allowed to change
! ::: dx          => cell size
! ::: problo      => phys loc of lower left corner of prob domain
! ::: time        => problem evolution time
! ::: level       => refinement level of this array
! ::: -----------------------------------------------------------

subroutine state_error(tag,tag_lo,tag_hi, &
                       state,state_lo,state_hi, &
                       set,clear,&
                       lo,hi,&
                       dx,problo,time,nc,maxgrox,maxgroy, tagfrac) bind(C, name="state_error")

  use amrex_fort_module, only : amrex_real
  use amrex_paralleldescriptor_module
  ! use amrex_parallel_module
  use amrex_amr_module
  implicit none
  
  integer          :: lo(3),hi(3),nc
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  real(amrex_real) :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),0:nc-1)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  real(amrex_real) :: problo(3), dx(3), time, maxgrox, maxgroy, tagfrac
  integer          :: set,clear

  integer          :: i, j, k, rank
  real(amrex_real) :: gradro, gro_x, gro_y, xcm, ycm, mid(2)
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  rank = amrex_pd_myproc()
  ! print*,"rank= ",rank,", fortran maxgradp= ",maxgradp,"lo= ",lo,", hi= ",hi
  ! gradpmax = abs(maxgradp)
  ! print*,"maxgradp= ",gradpmax
  ! print*, "maxgradro= ",maxgradro

  ! Locate the center of the square wave (find distance traveled by wave at time t)
  ! mid gives the x and y coordinates of the centre of the wave at t=0
  ! mid(1) = 0.5_amrex_real*(problo(1) + probhi(1))
  ! mid(2) = 0.5_amrex_real*(problo(2) + probhi(2))
  ! xcm = mid(1) + u*time
  ! ycm = mid(2) + v*time
  ! Tag on regions of high phi
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
            gro_x = 0.5_amrex_real*(state(i+1,j,k,ro) - state(i-1,j,k,ro))/dx(1)
            gro_y = 0.5_amrex_real*(state(i,j+1,k,ro) - state(i,j-1,k,ro))/dx(2)
            ! gradro = sqrt(gro_x**2 + gro_y**2)
            if(abs(gro_x) > tagfrac*maxgrox .or. &
            &  abs(gro_y) > tagfrac*maxgroy) then
            ! if (abs(gradro) .ge. 0.1_amrex_real*maxgradro) then 
            ! if (abs(gradro) > 0.05_amrex_real*maxgradro) then 
              tag(i,j,k) = set
            endif
        enddo
     enddo
  enddo

end subroutine state_error

