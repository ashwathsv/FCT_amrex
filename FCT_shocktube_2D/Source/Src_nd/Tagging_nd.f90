
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
                       dx,problo,time,phierr,nc,domdir) bind(C, name="state_error")

  use amrex_fort_module, only : amrex_real
  implicit none
  
  integer          :: lo(3),hi(3),nc,domdir
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  real(amrex_real) :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),0:nc-1)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  real(amrex_real) :: problo(3),dx(3),time,phierr
  integer          :: set,clear

  integer          :: i, j, k
  real(amrex_real) :: maxgradp, gradp

  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  ! We are tagging based on absolute value of pressure gradient
  ! First get the maximum pressure gradient
  maxgradp = 0.0_amrex_real
  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        if (domdir == 1) then
          gradp = 0.5_amrex_real*(state(i+1,j,k,pre) - state(i-1,j,k,pre))/dx(1)
        else
          gradp = 0.5_amrex_real*(state(i,j+1,k,pre) - state(i,j-1,k,pre))/dx(2)
        endif
        if (abs(gradp) > abs(maxgradp)) then
          maxgradp = gradp
        endif
      enddo
    enddo
  enddo
  print*,"maxgradp= ",maxgradp

  ! Tag on regions of high phi
  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
          if(domdir == 1) then
            gradp = 0.5_amrex_real*(state(i+1,j,k,pre) - state(i-1,j,k,pre))/dx(1)
            if (abs(gradp) .ge. 0.5_amrex_real*abs(maxgradp) .and. &
            &  i > lo(1)+4 .and. i < hi(1)-4) then
              tag(i,j,k) = set
            endif
         else
            gradp = 0.5_amrex_real*(state(i,j+1,k,pre) - state(i,j-1,k,pre))/dx(2)
            if (abs(gradp) .ge. 0.5_amrex_real*abs(maxgradp) .and. &
            &  j > lo(2)+4 .and. j < hi(2)-4) then
              tag(i,j,k) = set
            endif
         endif
        enddo
     enddo
  enddo

end subroutine state_error

