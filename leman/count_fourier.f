program main

  implicit none

  integer, parameter :: fn = 2
  integer mmin,mm,nmin,nn,n,m
  integer, dimension(:,:), allocatable :: mnTable

  character(253) s, stmp

  open(fn, file='AFtable.txt')
  read(fn, *)
  read(fn, *)
  read(fn, *)
  read(fn, '(1x,4i10)') mmin,mm,nmin,nn
  read(fn, *)
  read(fn, *)

  allocate(mnTable(1:mm,1:nn))

  write(stmp,'(i4)') mm
  s = '(5x,' // trim(stmp) // 'i2,i5)'
  s = trim(s)

  do n = 1,nn
    read(fn,s) (mnTable(m,n), m=1,mm)
  enddo

  close(fn)

  print *, 'Nb of AF modes:', sum(mnTable)

end program
