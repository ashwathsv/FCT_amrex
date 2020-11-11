
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
                       dx,problo,time,comp_lo,comp_hi,&
                       maxgradpx,maxgradpy,tagfrac,rad_bw, &
                       lev, max_level, lev_allow) bind(C, name="state_error")

  use amrex_fort_module, only : amrex_real
  ! use amrex_paralleldescriptor_module
  ! use amrex_parallel_module
  ! use amrex_amr_module
  implicit none
  
  integer          :: lo(3),hi(3),comp_lo,comp_hi,lev, max_level, lev_allow
  integer          :: state_lo(3),state_hi(3)
  integer          :: tag_lo(3),tag_hi(3)
  real(amrex_real) :: state(state_lo(1):state_hi(1), &
                            state_lo(2):state_hi(2), &
                            state_lo(3):state_hi(3),comp_lo:comp_hi)
  integer          :: tag(tag_lo(1):tag_hi(1),tag_lo(2):tag_hi(2),tag_lo(3):tag_hi(3))
  real(amrex_real) :: problo(3), dx(3), time, maxgradpx, maxgradpy, tagfrac, rad_bw
  integer          :: set,clear

  integer          :: i, j, k, rank
  real(amrex_real) :: dro, xcm, ycm, mid(2), gradpx, gradpy, dist
  integer, parameter :: ro = 0, rou = 1, rov = 2, roE = 3, pre = 4

  ! rank = amrex_pd_myproc()
  ! print*,"rank= ",rank,", rofs= ",rofs,"lo= ",lo,", hi= ",hi
  ! Tag on regions of high phi
  ! dist = sqrt((0.5*dx(1)**2) + (0.5*dx(2)**2))
  ! if(dist > rad_bw .and. time == 0.0) then
  !   print*,"lev= ",lev,", time= ",time, ", dist= ", dist, ", rad_bw= ", rad_bw
  !   ! if(lev == 0) then
  !     do k = lo(3),hi(3)
  !       do j = lo(2),lo(2)
  !         do i = lo(1),lo(1)
  !           tag(i,j,k) = set
  !         enddo
  !       enddo
  !     enddo
  !   ! endif
  ! else

  if(lev < lev_allow) then
    do       k = lo(3), hi(3)
      do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
            gradpx = 0.5_amrex_real*(state(i+1,j,k,pre) - state(i-1,j,k,pre))/dx(1)
            gradpy = 0.5_amrex_real*(state(i,j+1,k,pre) - state(i,j-1,k,pre))/dx(2)
            if(abs(gradpx) > tagfrac*maxgradpx .or. &
            &  abs(gradpy) > tagfrac*maxgradpy) then
              tag(i,j,k) = set
            endif
        enddo
      enddo
    enddo
  endif
  ! endif

end subroutine state_error

