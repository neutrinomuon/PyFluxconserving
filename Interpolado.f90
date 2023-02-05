! ###########################################################################
!     RESUME : For a given array of values 'xx_value' for the               !
!              abscissa, then return the ordinate array values of           !
!              'yy_value' based on a linear interpolation within a          !
!              table of pair values [xold_vec, yold_vec].                   !
!                                                                           !
!     Input           arguments = 5                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 2                                         !
!     Total number of arguments = 9                                         !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              03) xold_vec  -> Old x vector (abcissas)                     !
!              04) yold_vec  -> Old y vector (ordenadas)                    !
!              05) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              06) Is_Index  -> Search: Linear: [1] & Binary [0]            !
!              07) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points         !
!              02) IsKeepOn -> Flag, if == 0 then there's a problem         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Interpolado( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,IsKeepOn,Is_Index,verbosity )
    use ModDataType
    implicit none

    ! Variable declarations
    integer  (kind=IB), intent(in) :: nold_vec,nxyvalue
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB) :: i,j,k,ii,jj,indices_,IsShowOn
    integer  (kind=IB), optional :: Is_Index,verbosity
    
    real     (kind=RP), dimension(0:nold_vec-1), intent(in)  :: xold_vec,yold_vec
    real     (kind=RP), dimension(0:nxyvalue-1), intent(in)  :: xx_value
    real     (kind=RP), dimension(0:nxyvalue-1), intent(out) :: yy_value
    
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: nxyvalue=shape(yy_value,0)

    !f2py real     (kind=RP), intent(in) :: xold_vec
    !f2py real     (kind=RP), intent(in) :: yold_vec
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: nold_vec=shape(yold_vec,0)

    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), optional :: Is_Index=0
    !f2py integer  (kind=IB), optional :: verbosity=0

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
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[Interpolado]'
    end if

    if ( nold_vec < 1 ) then
       if ( IsShowOn == 1_IB ) then
          write (*,'(4x,a)') '[PROBLEM_INT] @@@@@@@@@@@@@@@@@@@@@@@@'
          write (*,'(4x,a)') '[Interpolado] nold_vec < 1 dimension @'
       end if
       yy_value = -999.0_RP    
       IsKeepOn = 0_IB
       return
    end if
    
    jj = 0
    do i=0,nxyvalue-1

       if ( indices_ == 1_IB ) then
          ! Linear Search
          do j=0,nold_vec-2
             if ( xx_value(i) >= xold_vec(j) .and.                          &
                  xx_value(i) <= xold_vec(j+1) ) then
                jj = j
                exit
             end if
          end do
       else
          ! Binary Search
          ii = 0
          if ( xx_value(i) > xold_vec(0) .and. xx_value(i) < xold_vec(nold_vec-1)  ) then
             ii = 0
             jj = nold_vec
             do while (jj > ii+1)
                k = (ii+jj)/2
                if (xx_value(i) < xold_vec(k)) then
                   jj=k
                else
                   ii=k
                end if
             end do
             jj = ii
          else if ( xx_value(i) <= xold_vec(0) ) then
             jj = 0
          else if ( xx_value(i) >= xold_vec(nold_vec-1) ) then
             jj = nold_vec-2
          end if
       end if

       if ( jj > nold_vec-2 ) then
          jj = nold_vec-2
       end  if
       if ( jj < 0 ) then
          jj = 0
       end  if
       
       yy_value(i) = yold_vec(jj)+(yold_vec(jj+1)-yold_vec(jj))             &
                   * (xx_value(i)-xold_vec(jj))/(xold_vec(jj+1)-xold_vec(jj))
    end do

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[Interpolado]'
    end if

    return
END SUBROUTINE Interpolado
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_Interpolado( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_Interpolado
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 002                                                          !
!
! 1) Interpolado
! 2) author_Interpolado
