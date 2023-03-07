!-------------------------------------------
subroutine CPUTime
  integer c,t,ic_max
  real, save :: stime
  integer, save :: icall
  real eltime
  
  
  call CPU_Time(eltime) !System_Clock(c,t,ic_max)

  if (icall == 0) then  ! Initial call
    stime = eltime
    call CPU_Time(eltime) !System_Clock(c,t,ic_max)
    icall = 1
  endif
  
  write(*,'(/a,1pe14.6,a)') 'CPU TIME USED SO FAR =',eltime-stime,' SECS'
  
end subroutine
!-------------------------------------------
