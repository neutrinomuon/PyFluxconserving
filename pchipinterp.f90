program test
  use ModDataType
  use pchip_module

  integer  (kind=IB), parameter :: NlAmax=20000
  real     (kind=RP), dimension(NlAmax) :: l,f,d,e,d_alt
  real     (kind=RP), dimension(:), allocatable :: l_new, f_new, f_new_spline

  integer  (kind=IB) :: m, findloc_min, findloc_max, Nl_max_new, n

  real     (kind=RP) :: xx

  integer  (kind=IB), dimension(2) :: next  !! integer array indicating number of extrapolation points:
  
  interface
     function deriv3(xx, xi, yi, ni, m)
       use ModDataType
       real     (kind=RP) :: deriv3, xx
       integer  (kind=IB) :: ni, m
       real     (kind=RP) :: xi(ni), yi(ni)
     end function deriv3
  end interface

  l     = -999.0
  e     = -999.0
  f     = -999.0
  d_alt = -999.0
  
  open  (unit=1,file='file.txt',status='old')
  do i=1,NlAmax
     read  (1,*,iostat=n,END=20) l(i),e(i)
  end do
  close (1)

20 Nl_max = i - 1

  incfd = 1
  !l(1) = 0.
  !delta = 0.01
  !do i = 2,Nl_max
  !   l(i) = l(i-1) + delta
  !end do

  !e(1:Nl_max) = sin(l(1:Nl_max))**2
  
  f = 0.
  do i = 2,Nl_max
     f(i) = f(i-1) + e(i-1) 
     !write (*,*) l(i),e(i),f(i)
  end do

  m = 1
  do i = 1,Nl_max
     xx = l(i)
     d(i) = deriv3(xx, l, f, Nl_max, m)
     !write (*,*) l(i),e(i),f(i),d(i)
  end do

  !stop
  
  call DPCHIM( Nl_max, l(1:Nl_max), f(1:Nl_max), d_alt(1:Nl_max), incfd, ierr )
  !write (*,*) ierr
   
  do i = 1,Nl_max
     !write (*,*) l(i),e(i),f(i),d(i),d_alt(i)
  end do

  Nl_max_new = Nl_max
  delta      = 10.
   
  allocate( l_new(Nl_max_new) )
  allocate( f_new(Nl_max_new) )
  allocate( f_new_spline(Nl_max_new) )
   
  l_new(1) = l(1)
  do i = 2,Nl_max_new
      l_new(i) = l_new(i-1) + delta
  end do
   
  next = 0

  f_new = 0.
  f_new_spline = 0.
  
  do i=1,Nl_max_new
      findloc_min = minloc( ARRAY=l-l_new(i), DIM=1, MASK=l>=l_new(i) )
      findloc_max = maxloc( ARRAY=l-l_new(i), DIM=1, MASK=l<=l_new(i) )
      !write (*,*) findloc_min, findloc_max, l_new(i), l(1), l(Nl_max)
   
      if ( (findloc_min > 0) .and. (findloc_max > 0) .and. (findloc_min /= findloc_max) ) then
         CALL dchfev ( l(findloc_min), l(findloc_max), f(findloc_min), f(findloc_max), &
                       d_alt(findloc_min), d_alt(findloc_max), 1, l_new(i), f_new(i), next, ierr)
      else
         f_new(i) = 0.
      end if
   
  end do
   
  ! Cubic SPLINE
  call SPLINE3DArr( l_new,f_new_spline,Nl_max_new,l,f,Nl_max,0 )
   
  do i = 1,Nl_max_new
      write (*,*) l_new(i),f_new(i),f_new_spline(i)
   end do
   
  deallocate( l_new )
  deallocate( f_new )
  deallocate( f_new_spline )
   
  stop
  ! -- 
  ! -- do i = 1,Nl_max
  ! --    if (i <= Nl_max/2) then
  ! --       write (*,*) l(i),e(i),f(i),l_new(i),f_new(i),d(i)
  ! --    else
  ! --       write (*,*) l(i),e(i),f(i),"0.000","0.000",d(i)
  ! --    end if
  ! -- end do

  
  
end program test
