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
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINECoeff( x,y,b,c,d,n,IsKeepOn,verbosity )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: verbosity
    
    integer  (kind=IB) :: i,j,IsShowOn

    real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
    real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d

    real     (kind=RP), dimension(0:n-2) :: h,M
    real     (kind=RP), dimension(0:n-1) :: l,u,z!,a

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

    IsKeepOn = 1_IB
    
    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB  ) then
       write (*,'(4x,a)') '[SPLINECoeff]'
    end if

    ! Part 1
    h(0:n-2) = x(1:n-1) - x(0:n-2)
    
    ! Part 2
    M(0) = 0
    M(1:n-2) = 3 * ( y(2:n-1) - y(1:n-2) ) / h(1:n-2) - 3 * ( y(1:n-2) - y(0:n-3) ) / h(0:n-3)
    
    ! Part 3
    l(0) = 1
    u(0) = 0
    z(0) = 0

    ! Part 4
    do i=1,n-2
       l(i) = 2 * ( x(i+1) - x(i-1) ) - h(i-1) * u(i-1)
       u(i) = h(i) / l(i)
       z(i) = ( M(i) - h(i-1) * z(i-1) ) / l(i)
    end do

    ! Part 5
    l(n-1) = 1
    z(n-1) = 0
    c(n-1) = 0

    ! Part 6
    do j=n-2,0,-1
       c(j) = z(j) - u(j) * c(j+1)
       b(j) = (y(j+1) - y(j)) / h(j) - h(j) * (c(j+1) + 2 * c(j)) / 3
       d(j) = (c(j+1) - c(j)) / (3 * h(j))
    end do
    
    ! Part 7 - Boundary Conditions    
    c(n-1) = 0
    b(n-1) = (y(n-1) - y(n-2)) / h(n-2) + h(n-2) * (c(n-2) + 2 * c(n-1)) / 3
    d(n-1) = (c(n-1) - c(n-2)) / (3 * h(n-2))

    if ( IsShowOn == 1_IB  ) then
       write (*,'(4x,a,10(e12.5))') '... b = ',(b(i),i=1,min(10,n))
       write (*,'(4x,a,10(e12.5))') '... c = ',(c(i),i=1,min(10,n))
       write (*,'(4x,a,10(e12.5))') '... d = ',(d(i),i=1,min(10,n))
       write (*,'(4x,a)') '[SPLINECoeff]'
    end if
    
    return
END SUBROUTINE SPLINECoeff
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
!     OUTPUT : 01) IfSPLINE3D -> Interpolated value at u                    !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION IfSPLINE3D( u,x,y,b,c,d,n,ilastval )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n,ilastval
    integer  (kind=IB) :: i,j,k,k1,k2,ilastnum,jlastnum
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
    !f2py real     (kind=RP), intent(in)  :: ilastval
    !f2py real     (kind=RP), intent(out)  :: IfSPLINE3D

    save ilastnum,jlastnum
    data ilastnum/0/
    data jlastnum/0/

    if ( ilastval <= 0_IB ) then
        i = 0
        j = n
     else
        i = ilastnum-1
        if ( u < x(jlastnum) ) then
           j = jlastnum+1
        else
           j = n
        end if
     end if
     !write(*,*) ilastnum,jlastnum!,!ilastval,i,j!,u
     
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    !i = 0
    if ( u > x(0) .AND. u < x(n-1)  ) then
       !i = 0
       !j = n
       !write (*,*) i,j

       !k = 0
       do while ( j > i+1 )

          !k1 = i + (j-i) / 3
          !k2 = j - (j-i) / 3

          !write (*,*) 'k1,k2 ',i,j,k1,k2,k
          
          !if ( u < x(k1) ) then
          !   j = k1
          !else
          !   if ( u < x(k2) ) then
          !      i = k1-1
          !      j = k2
          !   else
          !      i = k2
          !   end if
          !end if

          k = (i+j)/2
          if ( u < x(k) ) then
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
       j = 0
       !I_SPLINE3D = y(0)
       !return
    end if
    if ( u >= x(n-1) ) then
       i = n-1
       j = n-1
       !I_SPLINE3D = y(n-1)
       !return
    end if
! *** if u is ouside the x interval take a boundary value (left or right) ***    

    ilastnum = i
    jlastnum = j
    
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
    
! *** Evaluate spline interpolation *****************************************
    dx = u - x(i)
    IfSPLINE3D = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
    
    return
END FUNCTION IfSPLINE3D
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
!                                                                           !
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
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL (KIND=RP) FUNCTION SPLINE_func( u,x,y,n,IsKeepOn,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: IsShowOn,ilastval
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
       subroutine SPLINECoeff( x,y,b,c,d,n,IsKeepOn,verbosity )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d
       end subroutine SPLINECoeff
       real (kind=RP) function IfSPLINE3D( u,x,y,b,c,d,n,ilastval )
         use ModDataType
         integer  (kind=IB), intent(in) :: n,ilastval
         real     (kind=RP), intent(in), dimension(0:n-1) :: b,c,d
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(in) :: u
       end function IfSPLINE3D
    end interface

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if
    
! --- First compute the coefficients ----------------------------------------
    call SPLINECoeff( x,y,b,c,d,n,IsKeepOn,IsShowOn )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    ilastval = -999_IB
    SPLINE_Func = IfSPLINE3D( u,x,y,b,c,d,n,ilastval )
! --- Compute the output value ----------------------------------------------

    return
END FUNCTION SPLINE_Func
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
!                                                                           !
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
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINE3DAlt( u,o,m,x,y,n,IsKeepOn,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n,m
    integer  (kind=IB), intent(out) :: IsKeepOn

    integer  (kind=IB) :: i,IsShowOn,ilastval
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
       subroutine SPLINECoeff( x,y,b,c,d,n,IsKeepOn,verbosity )
         use ModDataType
         integer  (kind=IB), intent(in) :: n
         integer  (kind=IB), intent(out) :: IsKeepOn
         integer  (kind=IB), optional :: verbosity
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(out), dimension(0:n-1) :: b,c,d
       end subroutine SPLINECoeff
       real (kind=RP) function IfSPLINE3D( u,x,y,b,c,d,n,ilastval )
         use ModDataType
         integer  (kind=IB), intent(in) :: n,ilastval
         real     (kind=RP), intent(in), dimension(0:n-1) :: b,c,d
         real     (kind=RP), intent(in), dimension(0:n-1) :: x,y
         real     (kind=RP), intent(in) :: u
       end function IfSPLINE3D
    end interface

    IsKeepOn = 1_IB

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0
    end if
    
! --- First compute the coefficients ----------------------------------------
    call SPLINECoeff( x,y,b,c,d,n,IsKeepOn,IsShowOn )
! --- First compute the coefficients ----------------------------------------

! --- Compute the output value ----------------------------------------------
    if ( IsKeepOn == 1_IB ) then
       ilastval = -999_IB
       do i=0,m-1          
          o(i)=IfSPLINE3D( u(i),x,y,b,c,d,n,ilastval )

          if ( i==0 ) then
             ilastval = 1_IB
          end if

          if ( IsShowOn == 1_IB ) then
             if ( i == 0 ) then
                write (*,'(4x,a)') '[SPLINE3DAlt]'
             end if
             write (w,'(e15.8)') x(i)
             write (z,'(e15.8)') o(i)
             write (*,'(4x,a,a,a)') '... x(i): ',trim(adjustl(w))//' ==> y: ', &
                  trim(adjustl(z))
             if ( i == m-1 ) then
                write (*,'(4x,a)') '[SPLINE3DAlt]'
             end if
          end if
          
       end do
    else
       o = -999.0_RP
    end if
! --- Compute the output value ----------------------------------------------

    return    
END SUBROUTINE SPLINE3DAlt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_SPLINE3DAlt( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_SPLINE3DAlt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ###########################################################################
!     RESUME : Not working but leaving it here                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECURSIVE SUBROUTINE Ternary_Search( l, r, key, array, n, indice_l, indice_r )

    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB), intent(in) :: l,r
    integer  (kind=IB), intent(out) :: indice_l, indice_r

    integer  (kind=IB) :: mid1,mid2

    real     (kind=RP), intent(in), dimension(0:n-1) :: array
    real     (kind=RP), intent(in) :: key

    indice_l = l
    indice_r = r
    
    if ( r > l+1 ) then
       
       ! Find the mid1 and mid2
       mid1 = l + (r - l) / 3
       mid2 = r - (r - l) / 3
       
       if ( array(mid1) < key ) then
          indice_r = mid1
       else
          if ( array(mid2) < key ) then
             indice_l = mid1
             indice_r = mid2
          else
             indice_l = mid2
          end if
       end if
       
       call Ternary_Search( indice_l, indice_r, key, array, n, indice_l, indice_r )
       
    end if
    return
END SUBROUTINE Ternary_Search
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 006                                                          !
!
! 1) SPLINECoeff
! 2) IfSPLINE3D
! 3) SPLINE_func
! 4) SPLINE3DAlt
! 5) author_SPLINE3DAlt
! 6) Ternary_Search -> Not working
