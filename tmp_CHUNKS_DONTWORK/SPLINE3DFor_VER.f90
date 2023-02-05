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

    real     (kind=RP), intent(in), dimension(n) :: x,y
    real     (kind=RP), intent(out), dimension(n) :: b,c,d

    real     (kind=RP) :: h

    gap = n-1
! *** Check input values ****************************************************
    if ( n < 2 ) then
       return
    end if

    if ( n < 3 ) then
       b(1) = (y(2)-y(1))/(x(2)-x(1))                ! A linear interpolation
       c(1) = 0.0_RP
       d(1) = 0.0_RP
       b(2) = b(1)
       c(2) = 0.0_RP
       d(2) = 0.0_RP
       return
    end if

! *** Part 1: Setup values **************************************************
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
       d(i)   = x(i+1) - x(i)
       b(i)   = 2.0_RP*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i)   = c(i+1) - c(i)
    end do

! *** Part 2: Boundary conditions *******************************************
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0.0_RP
    c(n) = 0.0_RP

    if ( n /= 3 ) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

! *** Part 3: Forward elimination *******************************************

    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do

! *** 4: Back substitution **************************************************
    c(n) = c(n)/b(n)
    do j = 1, gap
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

! *** Part 5: Compute spline coefficients ***********************************
    b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
    do i = 1, gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.*c(i)
    end do
    c(n) = 3.0*c(n)
    d(n) = d(n-1)

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
    real     (kind=RP), intent(in), dimension(n) :: b,c,d
    real     (kind=RP), intent(in), dimension(n) :: x,y
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
    i = 1
    j = n+1
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
    I_SPLINE3D = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

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
         real     (kind=RP), intent(in), dimension(n) :: x,y
         real     (kind=RP), intent(out), dimension(n) :: b,c,d
       end subroutine SPLINE3DCoe
       real (kind=RP) function I_SPLINE3D( u,x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(n) :: b,c,d
         real     (kind=RP), intent(in), dimension(n) :: x,y
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
SUBROUTINE SPLINE3DArr( u,o,m,x,y,n,k )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n,m
    integer  (kind=IB) :: i
    real     (kind=RP), dimension(n) :: b,c,d
    real     (kind=RP), intent(in), dimension(n) :: x,y

    real     (kind=RP), intent(in), dimension(m) :: u
    real     (kind=RP), intent(out), dimension(m) :: o

    integer  (kind=IB), optional :: k
    
    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in)  :: u
    !f2py                     intent(hide), depend(u) :: m=shape(u,0)

    !f2py real     (kind=RP), intent(out)  :: o
    !f2py                     intent(hide), depend(o) :: m=shape(o,0)

    !f2py integer  (kind=IB), optional :: k=0
    
    interface
       subroutine SPLINE3DCoe( x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(n) :: x,y
         real     (kind=RP), intent(out), dimension(n) :: b,c,d
       end subroutine SPLINE3DCoe
       real (kind=RP) function I_SPLINE3D( u,x,y,b,c,d,n )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         real     (kind=RP), intent(in), dimension(n) :: b,c,d
         real     (kind=RP), intent(in), dimension(n) :: x,y
         real     (kind=RP), intent(in) :: u
       end function I_SPLINE3D
    end interface

! --- First compute the coefficients ----------------------------------------
    call SPLINE3DCoe( x,y,b,c,d,n )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    do i=1,m
       o(i)=I_SPLINE3D( u(i),x,y,b,c,d,n )
    end do
! --- Compute the output value ----------------------------------------------

END SUBROUTINE SPLINE3DArr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 001                                                          !
