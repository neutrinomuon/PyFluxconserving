! ###########################################################################
!     RESUME : Interpolation with spline 1D function. Evaluate the          !
!              coefficients b(i), c(i), and d(i), i=1,2,...,n for           !
!              cubic spline interpolation                                   !
!                                                                           !
!              s(x) = y(i) + b(i)*(x-x(i))**1 + c(i)*(x-x(i))**2            !
!                          + d(i)*(x-x(i))**3                               !
!                                                                           !
!              For x(i) <= x <= x(i+1)                                      !
!                                                                           !
!     INPUT  : 01) n -> # of elements in vector x and y                     !
!              02) x -> Old x vector (abcissas)                             !
!              03) y -> Old y vector (ordenadas)                            !
!                                                                           !
!     OUTPUT : 01) b -> Spline coefficient order 1                          !
!              02) c -> Spline coefficient order 2                          !
!              03) d -> Spline coefficient order 3                          !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINE3DCoe( x,y,b,c,d,n )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB) :: i,j,gap

    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
    real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d

    real     (kind=RP) :: h

    gap = n-2
! *** Check input values ****************************************************
    if ( n < 2 ) then
       return
    end if

    if ( n < 3 ) then
       b(0) = (y(1)-y(0))/(x(1)-x(0))                                       ! A linear interpolation
       c(0) = 0.0_RP
       d(0) = 0.0_RP
       b(1) = b(0)
       c(1) = 0.0_RP
       d(1) = 0.0_RP
       return
    end if

! *** Part 1: Setup values **************************************************
    d(0) = (x(1) - x(0))
    c(1) = (y(1) - y(0))/d(0)
    do i=1,gap
       d(i)   = x(i+1) - x(i)
       b(i)   = 2.0_RP*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i)   = c(i+1) - c(i)
    end do

! *** Part 2: Boundary conditions *******************************************
    b(0)   = -d(0)
    b(n-1) = -d(n-2)
    c(0)   = +0.0_RP
    c(n-1) = +0.0_RP

    if ( n /= 3 ) then
       c(0)   = c(2)/(x(3)-x(1)) - c(1)/(x(2)-x(0))
       c(n-1) = c(n-2)/(x(n-1)-x(n-3)) - c(n-3)/(x(n-2)-x(n-4))
       c(0)   = c(0)*d(0)**2/(x(3)-x(0))
       c(n-1) = -c(n-1)*d(n-2)**2/(x(n-1)-x(n-4))
    end if

! *** Part 3: Forward elimination *******************************************

    do i=1,n-1
       h    = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do

! *** 4: Back substitution **************************************************
    c(n-1) = c(n-1)/b(n-1)
    do j=0,gap
       i = n-j-1
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

! *** Part 5: Compute spline coefficients ***********************************
    b(n-1) = (y(n-1) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n-1))
    do i=0,gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.0_RP*c(i)
    end do
    c(n-1) = 3.0_RP*c(n-1)
    d(n-1) = d(n-2)

    !write (*,*) b,c,d
    
END SUBROUTINE SPLINE3DCoe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : This function evaluates the cubic spline interpolation       !
!              at point u, i.e.                                             !
!                                                                           !
!              ispline = y(i) + b(i)*(u-x(i))**1 + c(i)*(u-x(i))**2         !
!                             + d(i)*(u-x(i))**3                            !
!                                                                           !
!              Where x(i) <= u <= x(i+1)                                    !
!                                                                           !
!              The coefficients b(i), c(i), and d(i), i=1,2,...,n for       !
!              cubic spline interpolation is computed in the spline         !
!              subroutine.                                                  !
!                                                                           !
!              s(x) = y(i) + b(i)*(x-x(i))**1 + c(i)*(x-x(i))**2            !
!                          + d(i)*(x-x(i))**3                               !
!                                                                           !
!     INPUT  : 01) n -> # of elements in vector x and y                     !
!              02) x -> Old x vector (abcissas)                             !
!              03) y -> Old y vector (ordenadas)                            !
!              04) b -> Spline coefficient order 1                          !
!              05) c -> Spline coefficient order 2                          !
!              06) d -> Spline coefficient order 3                          !
!                                                                           !
!     OUTPUT : 01) u -> Interpolated y point                                !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION I_SPLINE3D( u,x,y,b,c,d,n )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB) :: i,j,k
    real     (kind=RP), intent(in), dimension(0:n-1) :: b,c,d
    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
    real     (kind=RP), intent(in) :: u
    real     (kind=RP) :: dx

! *** if u is ouside the x interval take a boundary value (left or right) ***
    !if ( u <= x(1) ) then
    !   I_SPLINE3D = y(1)
    !   return
    !end if
    !if ( u >= x(n) ) then
    !   I_SPLINE3D = y(n)
    !   return
    !end if

! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    i = 0
    j = n
    do while (j > i+1)
       k = (i+j)/2
       if(u < x(k)) then
          j=k
       else
          i=k
       end if
    end do

! *** Evaluate spline interpolation *****************************************
    dx = u - x(i)
    I_SPLINE3D = y(i) + dx*b(i) + dx**2*c(i) + dx**3*d(i)

!+ b(i)*(x-x(i))**1 + c(i)*(x-x(i))**2            !
!                          + d(i)*(x-x(i))**3        
    
    return
END FUNCTION I_SPLINE3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
REAL (KIND=RP) FUNCTION SPLINE3DFor( u,x,y,n )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    !integer  (kind=IB) :: i,j,k
    real     (kind=RP), dimension(n) :: b,c,d
    real     (kind=RP), intent(in), dimension(n) :: x,y
    real     (kind=RP), intent(in) :: u

    interface
       subroutine SPLINE3DCoe( x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d
       end subroutine SPLINE3DCoe
       real (kind=RP) function I_SPLINE3D( u,x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(0:n-1) :: b,c,d
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(in) :: u
       end function I_SPLINE3D
    end interface

! --- First compute the coefficients ----------------------------------------
    call SPLINE3DCoe( x,y,b,c,d,n )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    SPLINE3DFor = I_SPLINE3D( u,x,y,b,c,d,n )
! --- Compute the output value ----------------------------------------------

END FUNCTION SPLINE3DFor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
SUBROUTINE SPLINE3DArr( u,o,m,x,y,n,IsKeepOn,k )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n,m
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB) :: i,j,p,IsShowOn

    real     (kind=RP), allocatable, dimension(:) :: b,c,d
    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y

    real     (kind=RP), allocatable, dimension(:) :: res_x,res_y

    real     (kind=RP), intent(in), dimension(0:m-1) :: u
    real     (kind=RP), intent(out), dimension(0:m-1) :: o

    integer  (kind=IB), optional :: k
    
    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in)  :: u
    !f2py                     intent(hide), depend(u) :: m=shape(u,0)

    !f2py real     (kind=RP), intent(out)  :: o
    !f2py                     intent(hide), depend(o) :: m=shape(o,0)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), optional :: k=0
    
    interface
       subroutine SPLINE3DCoe( x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d
       end subroutine SPLINE3DCoe
       real (kind=RP) function I_SPLINE3D( u,x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(0:n-1) :: b,c,d
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(in) :: u
       end function I_SPLINE3D
    end interface

    IsKeepOn = 1_IB
    
    if ( present(k) ) then
        IsShowOn = k
    else
        IsShowOn = 0_IB
    end if
    
    if ( count(x(1:n-1)-x(0:n-2) < 0) > 0 ) then
       if ( IsShowOn == 1_IB ) then
          write (*,'(4x,a)') '... Houston, We have a problem'
       end if
       o = -999.0_RP
       IsKeepOn = 0_IB
       return
    end if

    if ( count(x(1:n-1)-x(0:n-2) == 0) > 0 ) then
       do i=0,n-2
          if ( x(i) == x(i+1) ) then
             !write (*,*) x(i),x(i+1)
             !write (*,*) y(i),y(i+1)

             if ( y(i) /= y(i+1) ) then
                if ( IsShowOn == 1_IB ) then
                   write (*,'(4x,a)') '... Houston, We have a problem'
                end if
                o = -999.0_RP
                IsKeepOn = 0_IB
                return
             end if
          end if
       end do

! *** Find unique values ****************************************************
       allocate( res_x(0:n-1) )
       allocate( res_y(0:n-1) )
       j = 0
       res_x    = 0.0_RP
       res_x(0) = x(0)
       res_y(0) = y(0)
       
       outer: do i=1,n-1
          p = 0
          inner: do while (p <= j)
             !write (*,*) i,p
             if ( res_x(p) == x(i) ) then
                !write (*,*) 'Olha que tem'
                cycle outer
             end if
       !      ! Found a match so start looking again
             p = p + 1
          end do inner
       !   ! No match found so add it to the output
          j = j + 1
          res_x(j) = x(i)
          res_y(j) = y(i)
       end do outer

       j = j + 1
       !write (*,*) j,n
       
       !deallocate( res_x )
       !deallocate( res_y )
! *** Find unique values ****************************************************

    else
       j = n
       allocate( res_x(0:n-1) )
       allocate( res_y(0:n-1) )
       res_x = x
       res_y = y
    end if
    
! --- First compute the coefficients ----------------------------------------
    allocate( b(0:j-1) )
    allocate( c(0:j-1) )
    allocate( d(0:j-1) )
    call SPLINE3DCoe( res_x(0:j-1),res_y(0:j-1),b,c,d,j )
! --- First compute the coefficients ----------------------------------------
    
! --- Compute the output value ----------------------------------------------
    do i=0,m-1
       o(i)=I_SPLINE3D( u(i),res_x(0:j-1),res_y(0:j-1),b,c,d,j )
    end do
! --- Compute the output value ----------------------------------------------

    deallocate( res_x )
    deallocate( res_y )

    deallocate( b )
    deallocate( c )
    deallocate( d )

    return
END SUBROUTINE SPLINE3DArr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 001                                                          !
