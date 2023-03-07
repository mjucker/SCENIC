      function sdot (dim,v1,d1,v2,d2)                                      12
C...Translated by  fopp     4.02E36 17:59:18  12/14/93
C...Switches: -eadejlpuvx18 -dchimr07 -e1 -gi
      integer j1
      integer   dim,d1,d2                                                  13
      real :: v1(dim),v2(dim)                                              14
      real*8 d3
      d3 = 0                                                               16
*VDIR NODEP
      do j1 = 1, dim                                                       16
         d3 = d3 + v1(1+(j1-1)*d1)*v2(1+(j1-1)*d2)                         16
      end do                                                               16
      sdot = d3                                                            16
      return                                                               17
      end                                                                  17
