! ###########################################################################
!     RESUME : For a given array of values 'xx_value' for the               !
!              abscissa, then return the ordinate array values of           !
!              'yy_value' based on a linear interpolation within a          !
!              table of pair values [xold_vec, yold_vec]. This              !
!              subroutine assumes that the values in the array              !
!              xold_vec increases monotonically with yold_vec.              !
!                                                                           !
!     Input           arguments = 5                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 1                                         !
!     Total number of arguments = 8                                         !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              03) xold_vec  -> Old x vector (abcissas)                     !
!              04) yold_vec  -> Old y vector (ordenadas)                    !
!              05) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              06) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points         !
!              02) IsKeepOn -> Flag, if == 0 then there's a problem         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LINdexerpol( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,IsKeepOn,verbosity )

    use ModDataType

    implicit none
    integer  (kind=IB), intent(in) :: nold_vec,nxyvalue
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: in_loop2,N_indice,IsShowOn

    real     (kind=RP), dimension(0:nold_vec-1), intent(in) :: xold_vec,yold_vec
    real     (kind=RP), dimension(0:nxyvalue-1), intent(in) :: xx_value
    real     (kind=RP), dimension(0:nxyvalue-1), intent(out) :: yy_value

    character (len=CH) :: W1aux,W2aux

    real     (kind=RP), dimension(0:nold_vec-1) :: a_vec,b_vec

    intrinsic adjustl, maxloc, trim
    
    !f2py real     (kind=RP), intent(in) :: xold_vec
    !f2py real     (kind=RP), intent(in) :: yold_vec
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: nold_vec=shape(yold_vec,0)

    !f2py real     (kind=RP), intent(in)  :: xx_value 
    !f2py real     (kind=RP), intent(out) :: yy_value 
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: nxyvalue=shape(yy_value,0)

    !f2py integer  (kind=IB), optional :: verbosity=0

    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
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
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[LINdexerpol]'
    end if
    
    a_vec(0:nold_vec-2)  = (yold_vec(1:nold_vec-1) - yold_vec(0:nold_vec-2))&
                         / (xold_vec(1:nold_vec-1) - xold_vec(0:nold_vec-2))
    b_vec(0:nold_vec-2)  =  yold_vec(0:nold_vec-2) - a_vec(0:nold_vec-2)    &
                         *  xold_vec(0:nold_vec-2)
    a_vec(nold_vec-1)    = (yold_vec(nold_vec-1)   - yold_vec(nold_vec-2))  &
                         / (xold_vec(nold_vec-1)   - xold_vec(nold_vec-2))
    b_vec(nold_vec-1)    =  yold_vec(nold_vec-1)   - a_vec(nold_vec-1)      &
                         *  xold_vec(nold_vec-1)
    
    N_indice = -999_IB
    do in_loop2=0,nxyvalue-1
       if ( xx_value(in_loop2) >= xold_vec(0) .and.                         &
            xx_value(in_loop2) <= xold_vec(nold_vec-1) ) then
          if ( N_indice < 0_IB ) then
             N_indice = maxloc( xold_vec, dim=1, mask=xold_vec<=xx_value(in_loop2) )
          else
             N_indice = N_indice                                            &
                      + maxloc( xold_vec(N_indice:nold_vec-1), dim=1,       &
                  mask=xold_vec(N_indice:nold_vec-1)<=xx_value(in_loop2) )
          end if
       else if ( xx_value(in_loop2) < xold_vec(0)  ) then
          N_indice = 1
       else if ( xx_value(in_loop2) > xold_vec(nold_vec-1)  ) then
          N_indice = nold_vec-1
       end if
       !write (*,*) N_indice
       !write (*,*) N_indice,xold_vec(N_indice),xx_value(in_loop2)
       !write (*,*) xx_value(in_loop2)
       
       yy_value(in_loop2) = a_vec(N_indice-1) * xx_value(in_loop2) + b_vec(N_indice-1)
    end do

    if ( IsShowOn == 1_IB ) then
       do in_loop2=1,nxyvalue
          write (W1aux,'(f15.4)') xx_value(in_loop2)
          write (W2aux,'(f15.4)') yy_value(in_loop2)
          write (*,'(4x,a,a,a,a)') '... xx_value(in_loop2): ',              &
               trim(adjustl(W1aux)),' ==> yy_value: ',trim(adjustl(W2aux))
       end do
    end if
       
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[LINdexerpol]'
    end if

    return
END SUBROUTINE LINdexerpol
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_LINdexerpol( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_LINdexerpol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Fri Sep 30 15:30:49 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!          PROGRAM GeneralTest
!            use ModDataType
!            implicit none
!            integer  (kind=IB) :: i,j,n
!            integer  (kind=IB), parameter :: Nl_max=80000
!            integer  (kind=IB) :: IsKeepOn,Int_Type,ilastval
!            real     (kind=RP), dimension(Nl_max) :: l,f,e,m
!            real     (kind=RP) :: int,lambda_i,lambda_f,ll1,ff1,ll2,ff2
!
!           END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 002                                                          !
!
! 1) LINdexerpol
! 2) author_LINdexerpol
