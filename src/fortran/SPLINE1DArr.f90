! ###########################################################################
!     RESUME : Spline1DArr subroutine is written in fortran 2003            !
!              standard. It takes an array of values x to                   !
!              interpolate from the arrays t and y.                         !
!                                                                           !
!     Input           arguments = 6                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 2                                         !
!     Total number of arguments = 10                                        !
!                                                                           !
!     INPUT  : 01) x        -> Interpolate into x array                     !
!              02) m        -> # of array elements in vector x and o        !
!              03) t        -> x vector (abcissas)                          !
!              04) y        -> y vector (ordenadas)                         !
!              05) n        -> # of array elements in vector t and y        !
!              06) e        -> epsilon machine dependent error              !
!              07) Is_Index -> Search: Linear: [1] & Binary [0]             !
!              08) k        -> Print & Check screen (OLD fortran way)       !
!                                                                           !
!     OUTPUT : 01) o        -> Interpolated y array                         !
!              02) IsKeepOn -> Flag, if == 0 then there's a problem         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!        Date: Sun Oct 15 10:00:52 WEST 2006                                !
!     Checked: Tue May  1 16:09:13 WEST 2012                                !
!              Fri Dec 28 14:55:10 WET  2012                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPLINE1DArr( x,o,m,t,y,n,e,IsKeepOn,Is_Index,verbosity )

    use ModDataType

    implicit none
    integer  (kind=IB), intent(in) :: n,m
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: Is_Index,verbosity
    integer  (kind=IB) :: i,j,k,IsShowOn
    integer  (kind=IB) :: ii,jj,indices_
    real     (kind=RP), intent(in), dimension(0:n-1)  :: t,y
    real     (kind=RP), intent(in), dimension(0:m-1)  :: x

    real     (kind=RP), intent(out), dimension(0:m-1) :: o
    
    real     (kind=RP), optional :: e

    real     (kind=RP) :: diff, epsilon
    character (len=CH) :: w,z,r

    !f2py real     (kind=RP), intent(in) :: t,y
    !f2py                     intent(hide), depend(t) :: n=shape(t,0)
    !f2py                     intent(hide), depend(y) :: n=shape(y,0)

    !f2py real     (kind=RP), intent(in) :: x
    !f2py                     intent(hide), depend(x) :: m=shape(x,0)

    !f2py real     (kind=RP), intent(out) :: o
    !f2py                     intent(hide), depend(o) :: m=shape(o,0)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    
    !f2py real     (kind=RP), intent(in), optional :: e=1.0e-8
    !f2py                     intent(in), optional :: Is_Index=0
    !f2py                     intent(in), optional :: verbosity=0
    
    intrinsic adjustl, count, present, trim

    if ( present(verbosity) ) then
       IsShowOn = verbosity
    else
       IsShowOn = 0_IB
    end if

    if ( present(Is_Index) ) then
       indices_ = Is_Index
    else
       indices_ = 0_IB
    end if

    IsKeepOn = 1_IB
    if ( n < 1 ) then
       if ( IsShowOn == 1_IB ) then
          write (*,'(4x,a)') '[PROBLEM_INT] @@@@@@@@@@@@@@@@@@@@@@@@'
          write (*,'(4x,a)') '[SPLINE1DArr] n < 1 dimension @@@@@@@@'
       end if
       o        = -999.0_RP    
       IsKeepOn = 0_IB
       return
    end if
    
    if ( present(e) ) then
       if ( e > 0.0_RP ) then
          epsilon = e
       else
          epsilon = 1.0e-8_RP
       end if
    else
       epsilon = 1.0e-8_RP
    end if

    i = 0
    do j=0,m-1
    
       if ( IsShowOn == 1_IB ) then
          if ( count(t(0:n-1) >= x(j) - epsilon) == 0_IB ) then
             write (*,'(4x,a)') '[PROBLEMATIC] @@@@@@@@@@@@@@@@@@@@@@@@'
             write (*,'(4x,a)') '[SPLINE1DArr] Out of allowed range @@@'
             write (*,'(4x,a)') '[SPLINE1DArr] SPLINE1DArr set = 0.0 @@'
             write (w,'(e15.6)') x(j)
             write (r,'(e15.6)') t(0)
             write (z,'(e15.6)') t(n-1)
             write (*,'(4x,a,a6,a,a6,a,a6)') '[SPLINE1DArr] ',trim(adjustl(w)),' => ',trim(adjustl(r)),' -- ',trim(adjustl(z))
             ! Removed the zero values !
             !o(j) = 0.0_RP
          else if ( count(t(0:n-1) <= x(j) + epsilon) == 0_IB ) then
             write (*,'(4x,a)') '[PROBLEMATIC] @@@@@@@@@@@@@@@@@@@@@@@@'
             write (*,'(4x,a)') '[SPLINE1DArr] Out of allowed range @@@'
             write (*,'(4x,a)') '[SPLINE1DArr] SPLINE1DArr set = 0.0 @@'
             write (w,'(e15.6)') x(j)
             write (r,'(e15.6)') t(0)
             write (z,'(e15.6)') t(n-1)
             write (*,'(4x,a,a6,a,a6,a,a6)') '[SPLINE1DArr] ',trim(adjustl(w)),' => ',trim(adjustl(r)),' -- ',trim(adjustl(z))
             ! Removed the zero values !
             !o(j) = 0.0_RP
          end if
       end if

       if ( indices_ == 0_IB ) then
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
          !i = 0
          ii = 0
          if ( x(j) > t(0) .AND. x(j) < t(n-1)  ) then
             ii = 0
             jj = n
             do while (jj > ii+1)
                k = (ii+jj)/2
                if (x(j) < t(k)) then
                   jj=k
                else
                   ii=k
                end if
             end do
             i = ii
          else if ( x(j) <= t(0) ) then
             i = 0
          else if ( x(j) >= t(n-1) ) then
             i = n-2
          end if
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************

       else
! *** Normal Linear Search **************************************************
          do i=n-1,0,-1
             diff = x(j) - t(i)
             if (diff >= 0.0_RP) exit
          end do
! *** Normal Linear Search **************************************************
       end if
    
       if ( i > n-2 ) then
          i = n-2
       end if
       if ( i < 0 ) then
          i = 0
       end if

       !write (*,*) i
          
       o(j) = y(i) + (x(j) - t(i))*(y(i+1) - y(i))/(t(i+1) - t(i))
       !end if

       if ( IsShowOn == 1_IB ) then
          if ( j == 0) then
             write (*,'(4x,a)') '[SPLINE1DArr]'
          end if
          write (w,'(e15.8)') x(j)
          write (z,'(e15.8)') o(j)
          write (*,'(4x,a,a,a)') '... x(j): ',trim(adjustl(w))//' ==> y: ', &
               trim(adjustl(z))
          if ( j == m-1 ) then
             write (*,'(4x,a)') '[SPLINE1DArr]'
          end if
       end if

    end do

    return
END SUBROUTINE SPLINE1DArr
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_SPLINE1Darr( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_SPLINE1Darr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************
  
! *** Number : 002                                                          !
!
! 1) SPLINE1DArr
! 2) author_SPLINE1Darr
