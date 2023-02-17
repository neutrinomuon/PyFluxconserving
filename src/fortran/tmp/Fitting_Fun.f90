! ###########################################################################
!     RESUME : For a given array of values 'xx_value' for the               !
!              abscissa, then return the ordinate array values of           !
!              'yy_value' based on a polynomial interpolation using         !
!              LAPACK library and the old pair of values                    !
!              [xold_vec,yold_vec].                                         !
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
!       OBS. : Needs to be compiled with flag --lapack                      !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Fitting_Fun

  use ModDataType
  
contains

  FUNCTION polyfitting( vx,vy,nx,d,IsKeepOn )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in)                  :: nx
    integer  (kind=IB), intent(in)                  :: d
    integer  (kind=IB), intent(out)                 :: IsKeepOn
    
    real     (kind=RP), dimension(d+1)              :: polyfitting
    real     (kind=RP), dimension(nx), intent(in)   :: vx,vy
 
    real     (kind=RP), dimension(:,:), allocatable :: X
    real     (kind=RP), dimension(:,:), allocatable :: XT
    real     (kind=RP), dimension(:,:), allocatable :: XTX
 
    integer  (kind=IB) :: i,j!,yydegree
 
    integer  (kind=IB) :: n,lda,lwork
    integer  (kind=IB) :: info
    integer  (kind=IB), dimension(:), allocatable :: ipiv
    real     (kind=RP), dimension(:), allocatable :: work

    !f2py real     (kind=RP), intent(in) :: vx,vy
    !f2py                     intent(hide), depend(vx) :: nx=shape(vx,0)
    !f2py                     intent(hide), depend(vy) :: nx=shape(vy,0)

    !f2py integer  (kind=IB), intent(in) :: d
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn

    IsKeepOn = 1_IB

    !if ( nx < d ) then
    !   yydegree = nx
    !end if
    
    n     = d+1
    lda   = n
    lwork = n
 
    allocate( ipiv(n)     )
    allocate( work(lwork) )
    allocate(   XT(n ,nx) )
    allocate(    X(nx,n ) )
    allocate(  XTX(n ,n ) )
 
    ! Prepare the matrix
    do i=0,d!yydegree
       do j=1,nx
          X(j,i+1) = vx(j)**i
       end do
    end do
 
    XT  = transpose(X)
    XTX = matmul(XT,X)
 
    ! Calls to LAPACK subs DGETRF and DGETRI
    call DGETRF( n,n,XTX,lda,ipiv,info )
    !write (*,*) info
    !write (*,*) XT
    !write (*,*) XTX
    !if ( info /= 0_IB ) then
       !print *, "problem"
    !   IsKeepOn = 0_IB
    !   polyfitting = -999.0_RP
    !   return
    !end if
    call DGETRI( n,XTX,lda,ipiv,work,lwork,info )
    !if ( info /= 0_IB ) then
       !print *, "problem"
       !write (*,*) 'Problem'
    !   IsKeepOn = 0_IB
    !   polyfitting = -999.0_RP
    !   return
    !end if

    polyfitting = matmul( matmul( XTX,XT ),vy )
    !if ( nx < d ) then
    !   polyfitting(nx+1:d+1) = 0.0_RP
    !end if
    
    deallocate( ipiv )
    deallocate( work )
    deallocate( X    )
    deallocate( XT   )
    deallocate( XTX  )

    return
  END FUNCTION polyfitting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! ###########################################################################
!     RESUME : For a given array of values 'xx_value' for the               !
!              abscissa, then return the ordinate array values of           !
!              'yy_value' based on a polynomial interpolation using         !
!              LAPACK library and the old pair of values                    !
!              [xold_vec,yold_vec].                                         !
!                                                                           !
!     Input           arguments = 6                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 1                                         !
!     Total number of arguments = 9                                         !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              03) xold_vec  -> Old x vector (abcissas)                     !
!              04) yold_vec  -> Old y vector (ordenadas)                    !
!              05) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              06) yydegree  -> Polynomial degree (integer)                 !
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
SUBROUTINE poly_interp( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,yydegree,IsKeepOn,verbosity )
    use ModDataType

    implicit none

    integer  (kind=IB), intent(in) :: nold_vec,nxyvalue
    integer  (kind=IB), intent(in) :: yydegree
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB), intent(out) :: IsKeepOn
    integer  (kind=IB) :: i,IsShowOn

    real     (kind=RP), dimension(0:nold_vec-1), intent(in)  :: xold_vec,   &
                                                                yold_vec
    real     (kind=RP), dimension(0:nxyvalue-1), intent(in)  :: xx_value
    real     (kind=RP), dimension(0:nxyvalue-1), intent(out) :: yy_value

    real     (kind=RP), dimension(yydegree+1) :: a

    !f2py real     (kind=RP), intent(in) :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: nxyvalue=shape(yy_value,0)

    !f2py real     (kind=RP), intent(in) :: xold_vec,yold_vec
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: nold_vec=shape(yold_vec,0)

    !f2py integer  (kind=IB), intent(in) :: yydegree
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn

    !f2py                     intent(in), optional :: verbosity=0
    
! *** Verbosity mode ********************************************************
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if
! *** Verbosity mode ********************************************************

    IsKeepOn = 1_IB

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[poly_interp]'
    end if

    a = 0.0
    a(1:min(nold_vec+1,yydegree+1)) = polyfitting( xold_vec,yold_vec,nold_vec,min(nold_vec,yydegree),IsKeepOn )

    !write (*,*) min(nold_vec,yydegree)
    
    !if ( yydegree > nold_vec ) then
    !   a(nold_vec+1:yydegree+1) = 0.0
    !end if

    !a = polyfitting( xold_vec,yold_vec,nold_vec,yydegree,IsKeepOn )
    
    !write (*,*) a
    
    if ( IsKeepOn == 1_IB ) then
       yy_value = 0.0_RP
       do i=0,yydegree
          yy_value = yy_value + a(i+1) * xx_value**i
       end do
    else
       yy_value = -999.0_RP
    end if

    if ( IsShowOn == 1_IB ) then
       write (*,*) a(0:max(yydegree,20) )
       write (*,'(4x,a)') '[poly_interp]'
    end if

    return
END SUBROUTINE poly_interp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_Fitting_Fun( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_Fitting_Fun

  
END MODULE Fitting_Fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Tue Sep 27 18:38:40 AZOST 2011 +++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 003                                                          !
!
! 1) polyfitting
! 2) poly_interp
! 3) author_Fitting_Fun