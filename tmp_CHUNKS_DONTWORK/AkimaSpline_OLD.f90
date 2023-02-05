! ###########################################################################
!     RESUME : Akima spline fitting interpolation. based on the paper       !
!              "A New Method of Interpolation and Smooth Curve Fitting      !
!              Based on Local Procedures" (HIROSHI AKIMA).                  !
!                                                                           !
!     INPUT  : 01) Nlambdas -> Number of points                             !
!              02) New_xval -> Abscissa for interpolation                   !
!              02) Oxvector -> The arrays of data abscissas                 !
!              03) Oyvector -> The arrays of data ordinates                 !
!              04) ErrCheck -> In case of error return ZERO                 !
!                                                                           !
!     OUTPUT : 01) New_yval -> New ordinate point                           !
!                                                                           !
!     OBS.   : It is assumed that Oxvector is monotonically                 !
!              increasing. There should be left at least 2 points on        !
!              each border to avoid at maximum border effects.              !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  9 13:22:46 WEST 2012                                !
!              Sun Mar  3 14:51:53 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AkimaSpline( New_xval,New_yval,Nrlambda,Oxvector,Oyvector,       &
                        N_lambda,ErrCheck,verbosity )

    use ModDataType

    implicit none
    integer  (kind=IB), intent(out) :: ErrCheck
    integer  (kind=IB), intent(in) :: N_lambda,Nrlambda
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: Nvecsize,Nvec_aux,IsShowOn,i,j
    real     (kind=RP), dimension(N_lambda), intent(in) :: Oxvector,Oyvector
    real     (kind=RP), allocatable, dimension(:) :: O_auxvec,Ozvector

    real     (kind=RP), intent(in), dimension(Nrlambda) :: New_xval
    real     (kind=RP), intent(out), dimension(Nrlambda) :: New_yval

    real     (kind=RP) :: a,b,c
    character (len=CH) :: W1aux

    !f2py real     (kind=RP), intent(out) :: ErrCheck
    !f2py real     (kind=RP), intent(in) :: Oxvector,Oyvector
    !f2py                     intent(hide), depend(Oxvector) :: N_lambda=shape(Oxvector,0)
    !f2py                     intent(hide), depend(Oyvector) :: N_lambda=shape(Oyvector,0)
    !f2py                     intent(in), optional :: verbosity=0

    !f2py real     (kind=RP), intent(in)  :: New_xval
    !f2py real     (kind=RP), intent(out) :: New_yval 
    !f2py                     intent(hide), depend(xx_value) :: Nrlambda=shape(New_xval,0)
    !f2py                     intent(hide), depend(yy_value) :: Nrlambda=shape(New_yval,0)

    intrinsic adjustl, present, trim
    
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    Nvecsize = N_lambda!size(Oxvector)
    Nvec_aux = Nvecsize+3
    allocate( O_auxvec(Nvec_aux) )
    allocate( Ozvector(Nvecsize) )

    O_auxvec = 0.0_RP
    Ozvector = 0.0_RP
    ErrCheck = 1_IB

    do j=1,Nrlambda

! *** Very special case where New_xval = 0 **********************************
       !if ( New_xval(j) == 0.0_RP ) then
       !   if ( IsShowOn == 1_IB ) then
       !      write (*,'(4x,a)') '[PROBLEM_BAS] @@@@@@@@@@@@@@@@@@@@@@@@'
       !      write (*,'(4x,a)') '[AkimaInterp] New_yval = 0.0 @@@@@@@@@'
       !      write (*,'(4x,a)') '[AkimaInterp] AkimaInterp set = 0.0 @@'
       !   end if
       !   New_yval(j) = 0.0_RP
       !   !ErrCheck = 0_IB
       !   !return
       !end if
! *** Very special case where New_xval = 0 **********************************

! *** Check if interpolation point is correct *******************************
       if ( New_xval(j) < Oxvector(1) ) then
          if ( IsShowOn == 1_IB ) then
             write (*,'(4x,a)') '[PROBLEM_BAS] @@@@@@@@@@@@@@@@@@@@@@@@'
             write (*,'(4x,a)') '[AkimaInterp] New_yval = 0.0 @@@@@@@@@'
             write (*,'(4x,a)') '[AkimaInterp] No Borders left @@@@@@@@'
             write (W1aux,'(f25.5)') New_xval
             write (*,'(4x,a,a)') '... New_xval: ',trim(adjustl(W1aux))
             write (W1aux,'(f25.5)') Oxvector(2)
             write (*,'(4x,a,a)') '... Oxvector: ',trim(adjustl(W1aux))
          end if
          New_yval(j) = Oyvector(1) + (Oyvector(2)-Oyvector(1))             &
                      * (New_xval(j)-Oxvector(1)) / (Oxvector(2)-Oxvector(1))
          ErrCheck = 0_IB
          !return
       end if
       if ( New_xval(j) > Oxvector(N_lambda) ) then
          if ( IsShowOn == 1_IB ) then
             write (*,'(4x,a)') '[PROBLEM_BAS] @@@@@@@@@@@@@@@@@@@@@@@@'
             write (*,'(4x,a)') '[AkimaInterp] New_yval = 0.0 @@@@@@@@@'
             write (*,'(4x,a)') '[AkimaInterp] No Borders right @@@@@@@'
             write (W1aux,'(f25.5)') New_xval
             write (*,'(4x,a,a)') '... New_xval: ',trim(adjustl(W1aux))
             write (W1aux,'(f25.5)') Oxvector(N_lambda-2)
             write (*,'(4x,a,a)') '... Oxvector: ',trim(adjustl(W1aux))
          end if
          New_yval(j) = Oyvector(N_lambda-1)                                &
                      + (Oyvector(N_lambda)-Oyvector(N_lambda-1))           &
                      * (New_xval(j)-Oxvector(N_lambda-1))                  &
                      / (Oxvector(N_lambda)-Oxvector(N_lambda-1))
          !New_yval(j) = 0.0_RP
          ErrCheck = 0_IB
          !return
       end if
! *** Check if interpolation point is correct *******************************

       if ( New_xval(j) >= Oxvector(1) .AND.                                &
            New_xval(j) <= Oxvector(N_lambda) ) then
          c = 2.0_RP * Oxvector(2) - Oxvector(3) ! Oxvector(1)
! *** Calculate Akima coefficients: a and b *********************************
!     RESUME : Shift i to i+2.                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i=2,N_lambda-1
             O_auxvec(i+2) = (Oyvector(i+1)-Oyvector(i))/(Oxvector(i+1)-Oxvector(i))
          end do

          O_auxvec(N_lambda+2) = 2.0_RP * O_auxvec(N_lambda+1) &
                               - O_auxvec(N_lambda)
          O_auxvec(N_lambda+3) = 2.0_RP * O_auxvec(N_lambda+2) &
                               - O_auxvec(N_lambda+1)
          O_auxvec(3)          = 2.0_RP * O_auxvec(4) - O_auxvec(5)
          O_auxvec(2)          = 2.0_RP * O_auxvec(3) - O_auxvec(4)

          i = 1
          a = 0.0_RP
          b = 0.0_RP
          do while ( i <= N_lambda .AND. a+b /=  0.0_RP )
             a = abs( O_auxvec(i+3)-O_auxvec(i+2) )
             b = abs( O_auxvec(i+1)-O_auxvec(i)   )
             
             if ( a+b /= 0.0_RP ) then
                Ozvector(i) = (a*O_auxvec(i+1)+b*O_auxvec(i+2))/(a+b)
             else
                Ozvector(i) = (O_auxvec(i+2)+O_auxvec(i+1))/2.0_RP
             end if
             
             i = i + 1
          end do

! *** Find table interval ***************************************************
          i = i + 1
          do while ( New_xval(j) > Oxvector(i) )
             i = i + 1
          end do
          i = i - 1
! *** Find table interval ***************************************************

          if ( i < 1 ) then
             i = 1
          end if
          if ( i > N_lambda ) then
             i = N_lambda-1
          end if
          
! *** Start interpolation scheme ********************************************
          b           = Oxvector(i+1) - Oxvector(i)
          a           = New_xval(j)   - Oxvector(i)
          New_yval(j) = Oyvector(i)  + Ozvector(i) * a + (3.0_RP*O_auxvec(i+2) &
                      - 2.0_RP * Ozvector(i) - Ozvector(i+1)) * a * a / b
          New_yval(j) = New_yval(j)  + (Ozvector(i)  + Ozvector(i+1)           &
                       - 2.0_RP * O_auxvec(i+2)) * a * a * a / (b * b)
! *** Start interpolation scheme ********************************************

          end if
          
    end do
 
    deallocate( O_auxvec )
    deallocate( Ozvector )
    return
END SUBROUTINE AkimaSpline
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_AkimaSpline( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_AkimaSpline
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 001                                                          !
