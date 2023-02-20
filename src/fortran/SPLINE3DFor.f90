! ###########################################################################
!     RESUME : Interpolation with cubic spline function. Evaluate the       !
!     coefficients b(i), c(i), and d(i), i=1,2,...,n for cubic spline       !
!     interpolation. In mathematics, a spline is a special function         !
!     defined piecewise by polynomials. In interpolating problems,          !
!     spline interpolation is often preferred to polynomial                 !
!     interpolation because it yields similar results, even when using      !
!     low degree polynomials, while avoiding Runge's phenomenon for         !
!     higher degrees.                                                       !
!                                                                           !
!              s(x) = y(i) + b(i)*(x-x(i))**1 + c(i)*(x-x(i))**2            !
!                          + d(i)*(x-x(i))**3                               !
!                                                                           !
!              For x(i) <= x <= x(i+1), strictly valid otherwise            !
!              extrapolate.                                                 !
!                                                                           !
!     Input           arguments = 3                                         !
!     Output          arguments = 4                                         !
!     Optional        arguments = 1                                         !
!     Total number of arguments = 8                                         !
!                                                                           !
!     INPUT  : 01) x         -> Old x vector (abcissas)                     !
!              02) y         -> Old y vector (ordenadas)                    !
!              03) n         -> # of elements in vector x and y             !
!              04) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) b         -> Spline coefficient order 1                  !
!              02) c         -> Spline coefficient order 2                  !
!              03) d         -> Spline coefficient order 3                  !
!              04) IsKeepOn  -> Flag, if == 0 then there's a problem        !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     References :                                                          !
!     1) Forsythe, G.E. (1977) Computer Methods For Mathematical            !
!     Computations. Ed. Prentice-Hall, Inc.                                 !
!                                                                           !
!     2) Bartels et al. (1998)                                              !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINE3DCoe( x,y,b,c,d,n,IsKeepOn,verbosity )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: verbosity
    
    integer  (kind=IB) :: i,j,last,IsShowOn

    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
    real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d

    real     (kind=RP) :: h

    character (len=CH) :: W1aux
    
    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(out)  :: b,c,d
    !f2py                     intent(hide), depend(b) :: n=shape(b,0)
    !f2py                     intent(hide), depend(c) :: n=shape(c,0)
    !f2py                     intent(hide), depend(d) :: n=shape(d,0)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), optional :: verbosity=0

    intrinsic adjustl, min, trim
    
    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB  ) then
       write (*,'(4x,a)') '[SPLINE3DCoe]'
    end if
    
    last = n-2

! *** Check input values ****************************************************
    if ( n < 2 ) then
       IsKeepOn = 0_IB
       b        = -999.0_RP
       c        = -999.0_RP
       d        = -999.0_RP
       if ( IsShowOn == 1_IB  ) then
          write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
          write (W1aux,'(i15)') n
          write (*,'(4x,a,a)') '... n < 2 ==',trim(adjustl(W1aux))
          write (*,'(4x,a)')   '... Set spline coefficients == -999 @'
       end if
       return
    end if

    if ( n < 3 ) then
       b(0) = (y(1)-y(0))/(x(1)-x(0))                                       ! A simple linear interpolation if less than 3 points
       c(0) = 0.0_RP
       d(0) = 0.0_RP
       b(1) = + b(0)
       c(1) = 0.0_RP
       d(1) = 0.0_RP
       return
    end if

! *** Part 1: Setup initial values ******************************************
    d(0) = (x(1) - x(0))
    c(1) = (y(1) - y(0))/d(0)

    do i=1,last
       d(i)   =        (x(i+1) - x(i))
       b(i)   = 2.0_RP*(d(i-1) + d(i))
       c(i+1) =        (y(i+1) - y(i))/d(i)
       c(i)   =        (c(i+1) - c(i))
    end do

! *** Part 2: Boundary conditions *******************************************
    b(0)   = -d(0)
    b(n-1) = -d(n-2)
    c(0)   = +0.0_RP
    c(n-1) = +0.0_RP

    if ( n /= 3 ) then
       c(0)   = +c(2)/(x(3)-x(1)) - c(1)/(x(2)-x(0))
       c(n-1) = +c(n-2)/(x(n-1)-x(n-3)) - c(n-3)/(x(n-2)-x(n-4))
       c(0)   = +c(0)*d(0)**2/(x(3)-x(0))
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

    ! Indexes from 1 to n comparison
    ! j = 1   => i = n-1
    ! j = n-1 => i = n-(n-1) = 1
    ! Indexes from 0 to n-1 should be
    ! j = 0   => i = n-2
    ! j = n-2 => i = 0
    do j=0,last
       i    = n-(j+2)
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

! *** Part 5: Compute spline coefficients ***********************************
!     RESUME: 3rd degree polynomial                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    b(n-1) = (y(n-1) - y(last))/d(last) + d(last)*(c(last) + 2.0_RP*c(n-1)) !/ 3
    do i=0,last
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0_RP*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.0_RP*c(i)
    end do
    c(n-1) = 3.0_RP*c(n-1)
    d(n-1) = d(n-2)

    if ( IsShowOn == 1_IB  ) then
       write (*,'(4x,a,10(e12.5))') '... b = ',(b(i),i=1,min(10,n))
       write (*,'(4x,a,10(e12.5))') '... c = ',(c(i),i=1,min(10,n))
       write (*,'(4x,a,10(e12.5))') '... d = ',(d(i),i=1,min(10,n))
       write (*,'(4x,a)') '[SPLINE3DCoe]'
    end if

    return
END SUBROUTINE SPLINE3DCoe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : This function evaluates the cubic spline interpolation       !
!              at point u, i.e.                                             !
!                                                                           !
!              I_SPLINE3D = y(i) + b(i)*(u-x(i))**1 + c(i)*(u-x(i))**2      !
!                         + d(i)*(u-x(i))**3                                !
!                                                                           !
!              Where x(i) <= u <= x(i+1)                                    !
!                                                                           !
!              The coefficients b(i), c(i), and d(i), i=1,2,...,n for       !
!              cubic spline interpolation is computed in the spline         !
!              subroutine.                                                  !
!                                                                           !
!     Input           arguments = 7                                         !
!     Output          arguments = 1                                         !
!     Optional        arguments = 0                                         !
!     Total number of arguments = 8                                         !
!                                                                           !
!     INPUT  : 01) u          -> New x value (abcissa)                      !
!              02) x          -> Old x vector (abcissas)                    !
!              03) y          -> Old y vector (ordenadas)                   !
!              04) b          -> Spline coefficient order 1                 !
!              05) c          -> Spline coefficient order 2                 !
!              06) d          -> Spline coefficient order 3                 !
!              07) n          -> # of elements in vector x and y            !
!                                                                           !
!     OUTPUT : 01) I_SPLINE3D -> Interpolated value at u                    !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
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

    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in)  :: b,c,d
    !f2py                     intent(hide), depend(b) :: n=shape(b,0)
    !f2py                     intent(hide), depend(c) :: n=shape(c,0)
    !f2py                     intent(hide), depend(d) :: n=shape(d,0)
    
    !f2py real     (kind=RP), intent(in)  :: u
    !f2py real     (kind=RP), intent(out)  :: I_SPLINE3D
    
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    i = 0
    if ( u > x(0) .AND. u < x(n-1)  ) then
       i = 0
       j = n
       do while (j > i+1)
          k = (i+j)/2
          if (u < x(k)) then
             j=k
          else
             i=k
          end if
       end do

!       i = 1
!       j = n+1
!20     k = (i+j)/2
!       if ( u  < x(k) ) j = k
!       if ( u >= x(k) ) i = k
!       if ( j  > i+1  ) go to 20
       
    end if
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    
! *** if u is ouside the x interval take a boundary value (left or right) ***
    if ( u <= x(0) ) then
       i = 0
       !I_SPLINE3D = y(0)
       !return
    end if
    if ( u >= x(n-1) ) then
       i = n-1
       !I_SPLINE3D = y(n-1)
       !return
    end if
! *** if u is ouside the x interval take a boundary value (left or right) ***    
    
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    
! *** Evaluate spline interpolation *****************************************
    dx = u - x(i)
    I_SPLINE3D = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))

    return
END FUNCTION I_SPLINE3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : This function evaluates the cubic spline interpolation       !
!              at u (real value), i.e.                                      !
!                                                                           !
!              SPLINE3DFor = y(i) + b(i)*(u-x(i))**1 + c(i)*(u-x(i))**2     !
!                          + d(i)*(u-x(i))**3                               !
!                                                                           !
!              Where stricly valid from x(i) <= u <= x(i+1)                 !
!                                                                           !
!              The coefficients b(i), c(i), and d(i), i=1,2,...,n for       !
!              cubic spline interpolation is computed in the spline         !
!              subroutine.                                                  !
!                                                                           !                                                                           !
!     Input           arguments = 4                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 1                                         !
!     Total number of arguments = 7                                         !
!                                                                           !
!     INPUT  : 01) u           -> New x value (abcissa)                     !
!              02) x           -> Old x vector (abcissas)                   !
!              03) y           -> Old y vector (ordenadas)                  !
!              04) n           -> # of elements in vector x and y           !
!              05) verbosity   -> Print & Check screen                      !
!                                                                           !
!     OUTPUT : 01) SPLINE3DFor -> Interpolated value of y at point u        !
!              02) IsKeepOn    -> Flag, if == 0 then there's a problem      !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SPLINE3DFor( u,x,y,n,IsKeepOn,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn
    real     (kind=RP), dimension(0:n-1) :: b,c,d
    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
    real     (kind=RP), intent(in) :: u

    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in)  :: u
    !f2py real     (kind=RP), intent(out)  :: SPLINE3DFor

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), optional :: verbosity=0

    interface
       subroutine SPLINE3DCoe( x,y,b,c,d,n,IsKeepOn,verbosity )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
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

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if
    
! --- First compute the coefficients ----------------------------------------
    call SPLINE3DCoe( x,y,b,c,d,n,IsKeepOn,IsShowOn )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    SPLINE3DFor = I_SPLINE3D( u,x,y,b,c,d,n )
! --- Compute the output value ----------------------------------------------

    return
END FUNCTION SPLINE3DFor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : This function evaluates the cubic spline interpolation       !
!              for an array u (j going from 0 to m), i.e.                   !
!                                                                           !
!              SPLINE3DArr = y(i) + b(i)*(u(j)-x(i))**1                     !
!                                 + c(i)*(u(j)-x(i))**2                     !
!                                 + d(i)*(u(j)-x(i))**3                     !
!                                                                           !
!              Where stricly valid from x(i) <= u(j) <= x(i+1)              !
!                                                                           !
!              The coefficients b(i), c(i), and d(i), i=1,2,...,n for       !
!              cubic spline interpolation is computed in the spline         !
!              subroutine.                                                  !
!                                                                           !                                                                           !
!     Input           arguments = 5                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 1                                         !
!     Total number of arguments = 8                                         !
!                                                                           !
!     INPUT  : 01) u           -> New x value (abcissa)                     !
!              02) m           -> # of elements in vector u and o           !
!              03) x           -> Old x vector (abcissas)                   !
!              04) y           -> Old y vector (ordenadas)                  !
!              05) n           -> # of elements in vector x and y           !
!              06) verbosity   -> Print & Check screen                      !
!                                                                           !
!     OUTPUT : 01) o           -> Interpolated array of ys at points u      !
!              02) IsKeepOn    -> Flag, if == 0 then there's a problem      !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINE3DArr( u,o,m,x,y,n,IsKeepOn,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n,m
    integer  (kind=IB), intent(out) :: IsKeepOn

    integer  (kind=IB) :: i,IsShowOn
    real     (kind=RP), dimension(0:n-1) :: b,c,d
    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y

    real     (kind=RP), intent(in), dimension(0:m-1) :: u
    real     (kind=RP), intent(out), dimension(0:m-1) :: o

    integer  (kind=IB), optional :: verbosity

    character (len=CH) :: w,z

    !f2py real     (kind=RP), intent(in)  :: x,y
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in)  :: u
    !f2py                     intent(hide), depend(u) :: m=shape(u,0)

    !f2py real     (kind=RP), intent(out)  :: o
    !f2py                     intent(hide), depend(o) :: m=shape(o,0)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), optional :: verbosity=0

    intrinsic adjustl, trim
    
    interface
       subroutine SPLINE3DCoe( x,y,b,c,d,n,IsKeepOn,verbosity )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
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

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0
    end if
    
! --- First compute the coefficients ----------------------------------------
    call SPLINE3DCoe( x,y,b,c,d,n,IsKeepOn,IsShowOn )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    if ( IsKeepOn == 1_IB ) then
       do i=0,m-1
          o(i)=I_SPLINE3D( u(i),x,y,b,c,d,n )

          if ( IsShowOn == 1_IB ) then
             if ( i == 0 ) then
                write (*,'(4x,a)') '[SPLINE3DArr]'
             end if
             write (w,'(e15.8)') x(i)
             write (z,'(e15.8)') o(i)
             write (*,'(4x,a,a,a)') '... x(i): ',trim(adjustl(w))//' ==> y: ', &
                  trim(adjustl(z))
             if ( i == m-1 ) then
                write (*,'(4x,a)') '[SPLINE3DArr]'
             end if
          end if
          
       end do
    else
       o = -999.0_RP
    end if
! --- Compute the output value ----------------------------------------------

END SUBROUTINE SPLINE3DArr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_SPLINE3DFor( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_SPLINE3DFor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 005                                                          !
!
! 1) SPLINE3DCoe
! 2) I_SPLINE3D
! 3) SPLINE3DFor
! 4) SPLINE3DArr
! 5) author_SPLINE3DFor
