! ###########################################################################
!     RESUME : Clenshaw–Curtis quadrature and Fejér quadrature are          !
!     methods for numerical integration, or "quadrature", that are          !
!     based on an expansion of the integrand in terms of Chebyshev          !
!     polynomials. Equivalently, they employ a change of variables          !
!     x=cos(theta) and use a discrete cosine transform (DCT)                !
!     approximation for the cosine series. Besides having                   !
!     fast-converging accuracy comparable to Gaussian quadrature            !
!     rules, Clenshaw–Curtis quadrature naturally leads to nested           !
!     quadrature rules (where different accuracy orders share points),      !
!     which is important for both adaptive quadrature and                   !
!     multidimensional quadrature (cubature).                               !
!                                                                           !
!     Briefly, the function f(x) to be integrated is evaluated at the       !
!     N extrema or roots of a Chebyshev polynomial and these values         !
!     are used to construct a polynomial approximation for the              !
!     function.                                                             !
!                                                                           !
!     The interval is assumed to be from [-1,1].                            !
!                                                                           !
!     The abscissas are the cosines of equally spaced angles between        !
!     180 and 0 degrees, including the endpoints.                           !
!                                                                           !
!     x(i) = cos( (order-i) * pi / (order - 1) )                            !
!                                                                           !
!     for the order == 1, we have:                                          !
!                                                                           !
!     x(0) = 0.0_RP                                                         !
!                                                                           !
!     If the value of order is increased in a sensible way, then the        !
!     new set of abscissas will include the old ones.  One such             !
!     sequence would be order(k) = 2*k+1 for k=0, 1, 2, ...                 !
!                                                                           !
!     For interpolation using Lagrange polynomials, the                     !
!     Clenshaw-Curtis abscissas could be used as a better choice than       !
!     regularly spaced abscissas.                                           !
!                                                                           !
!     Input           arguments = 1                                         !
!     Output          arguments = 1                                         !
!     Optional        arguments = 0                                         !
!     Total number of arguments = 2                                         !
!                                                                           !
!     INPUT  : 01) n -> Number of elements in array                         !
!                                                                           !
!     OUTPUT : 01) x -> Output array                                        !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!  Reference:                                                               !
!                                                                           !
!    [1] Charles Clenshaw, Alan Curtis, A Method for Numerical              !
!        Integration on an Automatic Computer, Numerische Mathematik,       !
!        Volume 2, Number 1, December 1960, pages 197-205.                  !
!                                                                           !
!    [2] Philip Davis, Philip Rabinowitz, Methods of Numerical              !
!        Integration, Second Edition, Dover, 2007, ISBN: 0486453391,        !
!        LC: QA299.3.D28.                                                   !
!                                                                           !
!    [3] Joerg Waldvogel, Fast Construction of the Fejer and                !
!        Clenshaw-Curtis Quadrature Rules, BIT Numerical Mathematics,       !
!        Volume 43, Number 1, 2003, pages 1-18.                             !
!                                                                           !
!    [4] https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature   !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CCabscissas ( x,n )
    use ModDataType    
    implicit none

    real     (kind=RP), parameter :: pi = 3.141592653589793_RP

    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB) :: i
    
    real     (kind=RP), dimension(0:n-1) :: theta
    real     (kind=RP), dimension(0:n-1), intent(out) :: x
    
    !f2py real     (kind=RP), intent(out) :: x
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)

    intrinsic real, cos
  
    if ( n  < 1_IB ) then
       write (*,'(4x,a)') '[PROBLEM_INT] @@@@@@@@@@@@@@@@@@@@@@@@'
       write (*,'(4x,a)') '[CCabscissas] Clenshaw–Curtis problem@'
       write (*,'(4x,a)') '[CCabscissas] dimension =====> n < 1 @'
       return
    end if
    
    if ( n == 1_IB ) then
       x(0) = 0.0_RP
       return
    end if
    
    do i=0,n-1
       theta(i) = real( n-i, kind=RP ) * pi / real( n-1, kind=RP )
    end do
    
    x(0:n-1) = cos( theta(0:n-1) )

    return
END SUBROUTINE CCabscissas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : The same as the previos routine, but for an interval         !
!     [A,B].                                                                !
!                                                                           !
!     The interval is assumed to be from [A,B].                             !
!                                                                           !
!     The abscissas are the cosines of equally spaced angles between        !
!     180 and 0 degrees, including the endpoints.                           !
!                                                                           !
!     x(i) = cos( (order-i) * pi / (order - 1) )                            !
!                                                                           !
!     for the order == 1, we have:                                          !
!                                                                           !
!     x(0) = 0.0_RP                                                         !
!                                                                           !
!     If the value of order is increased in a sensible way, then the        !
!     new set of abscissas will include the old ones.  One such             !
!     sequence would be order(k) = 2*k+1 for k=0, 1, 2, ...                 !
!                                                                           !
!     For interpolation using Lagrange polynomials, the                     !
!     Clenshaw-Curtis abscissas could be used as a better choice than       !
!     regularly spaced abscissas.                                           !
!                                                                           !
!     Input           arguments = 3                                         !
!     Output          arguments = 1                                         !
!     Optional        arguments = 0                                         !
!     Total number of arguments = 4                                         !
!                                                                           !
!     INPUT  : 01) n -> # of elements in x                                  !
!              02) a -> Initial point of interval                           !
!              03) b -> Final point of interval                             !
!                                                                           !
!     OUTPUT : 01) x -> Output array                                        !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!  Reference:                                                               !
!                                                                           !
!    [1] Charles Clenshaw, Alan Curtis, A Method for Numerical              !
!        Integration on an Automatic Computer, Numerische Mathematik,       !
!        Volume 2, Number 1, December 1960, pages 197-205.                  !
!                                                                           !
!    [2] Philip Davis, Philip Rabinowitz, Methods of Numerical              !
!        Integration, Second Edition, Dover, 2007, ISBN: 0486453391,        !
!        LC: QA299.3.D28.                                                   !
!                                                                           !
!    [3] Joerg Waldvogel, Fast Construction of the Fejer and                !
!        Clenshaw-Curtis Quadrature Rules, BIT Numerical Mathematics,       !
!        Volume 43, Number 1, 2003, pages 1-18.                             !
!                                                                           !
!    [4] https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature   !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CC_absc__ab( x,n,a,b )
    use ModDataType
    implicit none

    real     (kind=RP), parameter :: pi = 3.141592653589793_RP
    
    integer  (kind=IB), intent(in) :: n
    integer  (kind=IB) i
    
    real     (kind=RP), dimension(0:n-1), intent(out) :: x
    real     (kind=RP), dimension(0:n-1) :: theta
    
    real     (kind=RP), intent(in) :: a
    real     (kind=RP), intent(in) :: b
    
    !f2py real     (kind=RP), intent(in)  :: a,b
    !f2py real     (kind=RP), intent(out) :: x
    !f2py                     intent(hide), depend(x) :: n=shape(x,0)
    
    intrinsic real, cos
    
    if ( n < 1_IB ) then
       write (*,'(4x,a)') '[PROBLEM_INT] @@@@@@@@@@@@@@@@@@@@@@@@'
       write (*,'(4x,a)') '[CCabscissas] Clenshaw–Curtis problem@'
       write (*,'(4x,a)') '[CCabscissas] dimension =====> n < 1 @'
       return
    end if
    
    if ( n == 1_IB ) then
       x(0) = 0.5_RP * ( b + a )
       return
    end if
    
    do i=0,n-1
       theta(i) = real ( n-i, kind=RP ) * pi / real ( n-1, kind=RP )
    end do
    
    x(0:n-1) = 0.5_RP * ( ( b + a ) + ( b - a ) * cos( theta(1:n) ) )
    
    return
END SUBROUTINE CC_absc__ab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! F1_ABSCISSAS computes Fejer type 1 abscissas.
!
!  Discussion:
!
!    The interval is [ -1, +1 ].
!
!    The abscissas are the cosines of equally spaced angles, which
!    are the midpoints of N equal intervals between 0 and PI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
SUBROUTINE f1_abscissas ( n, x )

  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F1_ABSCISSAS - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( 2 * n - 2 * i + 1, kind = 8 ) * pi &
             / real ( 2 * n,             kind = 8 )
  end do

  x(1:n) = cos ( theta(1:n) )

  return
END SUBROUTINE f1_abscissas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! F1_ABSCISSAS_AB computes Fejer type 1 abscissas for the interval [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
subroutine f1_abscissas_ab ( a, b, n, x )
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F1_ABSCISSAS_AB - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = 0.5D+00 * ( b + a )
    return
  end if

  do i = 1, n
    theta(i) = real ( 2 * n - 2 * i + 1, kind = 8 ) * pi &
             / real ( 2 * n,             kind = 8 )
  end do

  x(1:n) = 0.5D+00 * ( ( b + a ) + ( b - a ) * cos ( theta(1:n) ) )

  return
END subroutine f1_abscissas_ab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!! F2_ABSCISSAS computes Fejer Type 2 abscissas.
!
!  Discussion:
!
!    The interval is [-1,+1].
!
!    The abscissas are the cosines of equally spaced angles.
!    The angles are computed as N+2 equally spaced values between 0 and PI,
!    but with the first and last angle omitted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Walter Gautschi,
!    Numerical Quadrature in the Presence of a Singularity,
!    SIAM Journal on Numerical Analysis,
!    Volume 4, Number 3, 1967, pages 357-362.
!
!    Joerg Waldvogel,
!    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
!    BIT Numerical Mathematics,
!    Volume 43, Number 1, 2003, pages 1-18.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE f2_abscissas ( n, x )

  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F2_ABSCISSAS - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if 

  if ( n == 1 ) then
    x(1) = 0.0D+00
    return
  else if ( n == 2 ) then
    x(1) = -0.5D+00
    x(2) =  0.5D+00
    return
  end if

  do i = 1, n
    theta(i) = real ( n + 1 - i, kind = 8 ) * pi &
             / real ( n + 1,     kind = 8 )
  end do

  x(1:n) = cos ( theta(1:n) )

  return
END SUBROUTINE f2_abscissas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! F2_ABSCISSAS_AB computes Fejer Type 2 abscissas for the interval [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine f2_abscissas_ab ( a, b, n, x )
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) :: pi = 3.141592653589793D+00
  real ( kind = 8 ) theta(n)
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'F2_ABSCISSAS_AB - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  do i = 1, n
    theta(i) = real ( n + 1 - i, kind = 8 ) * pi &
             / real ( n + 1,     kind = 8 )
  end do

  x(1:n) = 0.5D+00 * ( ( b + a ) + ( b - a ) * cos ( theta(1:n) ) )

  return
END subroutine f2_abscissas_ab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : This routine assumes a space of mxyvalue dimensions for      !
!              the yold_vec and yy_value. It does a Lagrange polynomial     !
!              interpolation. For a given array of values 'xx_value'        !
!              for the abscissa, then return the ordinate array values      !
!              of 'yy_value' based on a linear interpolation within a       !
!              table of pair values [xold_vec, yold_vec]. This              !
!              subroutine assumes that the values in the array              !
!              xold_vec increases monotonically with yold_vec. Note         !
!              that the user may request extrapolation.                     !
!                                                                           !
!              For each mxyvalue dimension, a polynomial of degree          !
!              (nold_vec - 1) is generated for the                          !
!              interpolation. WARNING: for most cases, the polynomial       !
!              begins to strongly oscillate when the number of              !
!              datapoints nold_vec increases. Therefore, it is not          !
!              well behaved for higher datapoints. This happens even        !
!              if the original data is well behaved. Ideal values of        !
!              nold_vec are <= ten (10).                                    !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              02) xold_vec  -> Old x vector (abcissas)                     !
!              03) yold_vec  -> Old y vector (ordenadas)                    !
!              04) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              05) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points and     !
!                              mxyvalue dimensions.                         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE intlagrange ( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,      &
                         nold_vec,mxyvalue,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: nold_vec,mxyvalue,nxyvalue
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB) :: i,j,IsShowOn
    
    real     (kind=RP), dimension(nold_vec), intent(in) :: xold_vec
    real     (kind=RP), dimension(mxyvalue,nold_vec), intent(in)  :: yold_vec    

    real     (kind=RP), dimension(nxyvalue), intent(in) :: xx_value
    real     (kind=RP), dimension(mxyvalue,nxyvalue), intent(out) :: yy_value

    real     (kind=RP), dimension(nold_vec,nxyvalue) :: l_interp

    character (len=CH) :: W1aux,W2aux
    
    !f2py real     (kind=RP), intent(in)  :: xold_vec
    !f2py real     (kind=RP), intent(in)  :: yold_vec
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: mxyvalue=shape(yold_vec,0),nold_vec=shape(yold_vec,1)
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: mxyvalue=shape(yy_value,0),nxyvalue=shape(yy_value,1)
    
    !f2py integer  (kind=IB), optional :: verbosity=0

    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[intlagrange]'
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  *** Compute the nold_vec Lagrange polynomials associated with           !
!      xold_vec(1:nold_vec) for the interpolation points                   !
!      xx_value(1:nxyvalue).                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call lagrangeval (xx_value, l_interp, nxyvalue, xold_vec, nold_vec )

!  *** Multiply                                                            !
!      yold_vec(1:mxyvalue,1:nold_vec) * l_interp(1:nold_vec,1:nxyvalue)   !
!      to get the interpolated yy_value(1:mxyvalue,1:nxyvalue) values.     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    yy_value(1:mxyvalue,1:nxyvalue) =                                      &
  matmul ( yold_vec(1:mxyvalue,1:nold_vec),l_interp(1:nold_vec,1:nxyvalue) )

    do j=1,mxyvalue
       where( xx_value(1:nxyvalue) < xold_vec(1) )
          yy_value(j,1:nxyvalue) = yold_vec(j,1) 
       end where
       where( xx_value(1:nxyvalue) > xold_vec(nold_vec) )
          yy_value(j,1:nxyvalue) = yold_vec(j,nold_vec) 
       end where
    end do
 
    if ( IsShowOn == 1_IB ) then
       do j=1,nxyvalue
          write (W1aux,'(f15.4)') xx_value(j)
          write (W2aux,'(10(e12.4,4x)))') (yy_value(i,j),i=1,min(mxyvalue,10))
          write (*,'(4x,a,a,a,a)') '... xx_value(1,j): ',                  &
               trim(adjustl(W1aux)),' ==> yy_value: ',trim(adjustl(W2aux))
       end do
    end if
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[intlagrange]'
    end if
    
  return
END SUBROUTINE intlagrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : This routine assumes a space of mxyvalue dimensions for      !
!              the yold_vec and yy_value. It does a piecewise linear        !
!              interpolation. For a given array of values 'xx_value'        !
!              for the abscissa, then return the ordinate array values      !
!              of 'yy_value' based on a linear interpolation within a       !
!              table of pair values [xold_vec, yold_vec]. This              !
!              subroutine assumes that the values in the array              !
!              xold_vec increases monotonically with yold_vec.  Note        !
!              that the user may request extrapolation.                     !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              02) xold_vec  -> Old x vector (abcissas)                     !
!              03) yold_vec  -> Old y vector (ordenadas)                    !
!              04) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              05) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points and     !
!                              mxyvalue dimensions.                         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE interp__lin( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,mxyvalue,checking,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: nold_vec
    integer  (kind=IB), intent(in) :: mxyvalue
    integer  (kind=IB), intent(in) :: nxyvalue
    
    integer  (kind=IB), optional :: checking,verbosity
    
    integer  (kind=IB) :: interp,i,j,IsShowOn
    integer  (kind=IB) :: left
    integer  (kind=IB) :: right
    
    real     (kind=RP), dimension(nold_vec), intent(in) :: xold_vec
    real     (kind=RP), dimension(nxyvalue), intent(in) :: xx_value
    real     (kind=RP), dimension(mxyvalue,nold_vec), intent(in) :: yold_vec
    real     (kind=RP), dimension(mxyvalue,nxyvalue),intent(out) :: yy_value

    real     (kind=RP) :: x

    logical :: r8vec_ascends_strictly

    character (len=CH) :: W1aux,W2aux
    
    !f2py real     (kind=RP), intent(in)  :: xold_vec
    !f2py real     (kind=RP), intent(in)  :: yold_vec
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: mxyvalue=shape(yold_vec,0),nold_vec=shape(yold_vec,1)
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: mxyvalue=shape(yy_value,0),nxyvalue=shape(yy_value,1)
    
    !f2py integer  (kind=IB), optional :: checking=0
    !f2py integer  (kind=IB), optional :: verbosity=0
    
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp__lin]'
    end if
    
    if ( checking == 1 .and. .not. r8vec_ascends_strictly ( nold_vec, xold_vec ) ) then
       write (*,'(4x,a)') '[PROBLEM_BAS] @@@@@@@@@@@@@@@@@@@@@@@@'
       write (*,'(4x,a)') '[interp__lin] - Fatal error!'
       write (*,'(4x,a)') '... Independent variable array is not strictly increasing.'
       yy_value(1:mxyvalue,1:nxyvalue) = -999.0_RP
       return
    end if
    
    do interp=1,nxyvalue
       
       x = xx_value(interp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *** Find the interval [ xold_vec(left),xold_vec(right) ] encompassing, or !
!     near to, x.                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call r8vec_bracket ( nold_vec, xold_vec, x, left, right )

       yy_value(1:mxyvalue,interp) = ( ( xold_vec(right) - x              ) &
                                       * yold_vec(1:mxyvalue,left)          &
                                     + ( x               - xold_vec(left) ) &
                                       * yold_vec(1:mxyvalue,right) )       &
                                     / ( xold_vec(right) - xold_vec(left) )
    end do

    if ( IsShowOn == 1_IB ) then
       do j=1,nxyvalue
          write (W1aux,'(f15.4)') xx_value(j)
          write (W2aux,'(10(e12.4,4x)))') (yy_value(i,j),i=1,min(mxyvalue,10))
          write (*,'(4x,a,a,a,a)') '... xx_value(1,j): ',                   &
               trim(adjustl(W1aux)),' ==> yy_value: ',trim(adjustl(W2aux))
       end do
    end if
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp__lin]'
    end if

    return
END SUBROUTINE interp__lin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : This routine assumes a space of mxyvalue dimensions for      !
!              the yold_vec and yy_value. It does a piecewise linear        !
!              interpolation. For a given array of values 'xx_value'        !
!              for the abscissa, then return the ordinate array values      !
!              of 'yy_value' based on a linear interpolation within a       !
!              table of pair values [xold_vec, yold_vec]. This              !
!              subroutine assumes that the values in the array              !
!              xold_vec increases monotonically with yold_vec.  Note        !
!              that the user may request extrapolation.                     !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              02) xold_vec  -> Old x vector (abcissas)                     !
!              03) yold_vec  -> Old y vector (ordenadas)                    !
!              04) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              05) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points and     !
!                              mxyvalue dimensions.                         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE interp_1lin( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,checking,verbosity )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: nold_vec
    integer  (kind=IB), intent(in) :: nxyvalue
    
    integer  (kind=IB), optional :: checking,verbosity
    
    integer  (kind=IB) :: interp,i,j,IsShowOn
    integer  (kind=IB) :: left
    integer  (kind=IB) :: right
    
    real     (kind=RP), dimension(nold_vec), intent(in) :: xold_vec
    real     (kind=RP), dimension(nold_vec), intent(in) :: yold_vec
    real     (kind=RP), dimension(nxyvalue), intent(in) :: xx_value
    real     (kind=RP), dimension(nxyvalue),intent(out) :: yy_value

    real     (kind=RP) :: x

    logical :: r8vec_ascends_strictly

    character (len=CH) :: W1aux,W2aux
    
    !f2py real     (kind=RP), intent(in)  :: xold_vec
    !f2py real     (kind=RP), intent(in)  :: yold_vec
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hide), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hide), depend(yold_vec) :: nold_vec=shape(yold_vec,0)
    !f2py                     intent(hide), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hide), depend(yy_value) :: nxyvalue=shape(yy_value,0)
    
    !f2py integer  (kind=IB), optional :: checking=0
    !f2py integer  (kind=IB), optional :: verbosity=0
    
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp_1lin]'
    end if
    
    if ( checking == 1 .and. .not. r8vec_ascends_strictly( nold_vec, xold_vec ) ) then
       write (*,'(4x,a)') '[PROBLEM_BAS] @@@@@@@@@@@@@@@@@@@@@@@@'
       write (*,'(4x,a)') '[interp__lin] - Fatal error!'
       write (*,'(4x,a)') '... Independent variable array is not strictly increasing.'
       yy_value(1:nxyvalue) = -999.0_RP
       return
    end if
    
    do interp=1,nxyvalue
       
       x = xx_value(interp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! *** Find the interval [ xold_vec(left),xold_vec(right) ] encompassing, or !
!     near to, x.                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       call r8vec_bracket ( nold_vec, xold_vec, x, left, right )

       yy_value(interp) = ( ( xold_vec(right) - x              )            &
                        * yold_vec(left)                                    &
                        + ( x               - xold_vec(left) )              &
                        * yold_vec(right) )                                 &
                        / ( xold_vec(right) - xold_vec(left) )
    end do

    if ( IsShowOn == 1_IB ) then
       do j=1,nxyvalue
          write (W1aux,'(f15.4)') xx_value(j)
          write (W2aux,'(10(e12.4,4x)))') yy_value(j)
          write (*,'(4x,a,a,a,a)') '... xx_value(1,j): ',                   &
               trim(adjustl(W1aux)),' ==> yy_value: ',trim(adjustl(W2aux))
       end do
    end if
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp_1lin]'
    end if

    return
END SUBROUTINE interp_1lin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     RESUME : This routine assumes a space of mxyvalue dimensions for      !
!              the yold_vec and yy_value. It does a nearest neighbor        !
!              interpolation. For a given array of values 'xx_value'        !
!              for the abscissa, then return the ordinate array values      !
!              of 'yy_value' based on a linear interpolation within a       !
!              table of pair values [xold_vec, yold_vec]. This              !
!              subroutine assumes that the values in the array              !
!              xold_vec increases monotonically with yold_vec.  Note        !
!              that the user may request extrapolation.                     !
!                                                                           !
!     INPUT  : 01) xx_value  -> New interpolated x array with points        !
!              02) nxyvalue  -> # of elements in xx_value and yy_value      !
!              02) xold_vec  -> Old x vector (abcissas)                     !
!              03) yold_vec  -> Old y vector (ordenadas)                    !
!              04) nold_vec  -> # of elements in xold_vec and yold_vec      !
!              05) verbosity -> Print & Check screen                        !
!                                                                           !
!     OUTPUT : 01) yy_value -> New interpolated y array with points and     !
!                              mxyvalue dimensions.                         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     Written: Jean Michel Gomes                                            !
!     Checked: Wed May  2 10:00:52 WEST 2012                                !
!              Fri Dec 28 13:53:36 WET  2012                                !
!              Sun Mar 10 10:05:03 WET  2013                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE interp_near( xx_value,yy_value,nxyvalue,xold_vec,yold_vec,       &
                        nold_vec,mxyvalue,verbosity )
  
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: nold_vec,mxyvalue,nxyvalue
    integer  (kind=IB), optional :: verbosity
    
    integer  (kind=IB) :: i,j,IsShowOn
    integer  (kind=IB) ::  r8vec_sorted_nearest
    
    real     (kind=RP), dimension(nold_vec), intent(in)           :: xold_vec
    real     (kind=RP), dimension(mxyvalue,nold_vec), intent(in)  ::  yold_vec

    real     (kind=RP), dimension(nxyvalue), intent(in)           :: xx_value
    real     (kind=RP), dimension(mxyvalue,nxyvalue), intent(out) :: yy_value

    character (len=CH) :: W1aux,W2aux
    
    !f2py real     (kind=RP), intent(in)  :: xold_vec
    !f2py real     (kind=RP), intent(in)  :: yold_vec
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: yy_value
    !f2py                     intent(hiden), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hiden), depend(yold_vec) :: mxyvalue=shape(yold_vec,0),nold_vec=shape(yold_vec,1)
    !f2py                     intent(hiden), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hiden), depend(yy_value) :: mxyvalue=shape(yy_value,0),nxyvalue=shape(yy_value,1)

    !f2py                   , optional   :: verbosity=0
    
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if

    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp_near]'
    end if
    
! *** Find the nearest index and rearrange terms ****************************
    do j=1,nxyvalue
       i = r8vec_sorted_nearest( nold_vec, xold_vec, xx_value(j) )
       yy_value(1:mxyvalue,j) = yold_vec(1:mxyvalue,i)
    end do
! *** Find the nearest index and rearrange terms ****************************

    if ( IsShowOn == 1_IB ) then
       do j=1,nxyvalue
          write (W1aux,'(f15.4)') xx_value(j)
          write (W2aux,'(10(e12.4,4x)))') (yy_value(i,j),i=1,min(mxyvalue,10))
          write (*,'(4x,a,a,a,a)') '... xx_value(1,j): ',                   &
               trim(adjustl(W1aux)),' ==> yy_value: ',trim(adjustl(W2aux))
       end do
    end if
    
    if ( IsShowOn == 1_IB ) then
       write (*,'(4x,a)') '[interp_near]'
    end if

    return
END SUBROUTINE interp_near
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LAGRANGE_VALUE evaluates the Lagrange polynomials.
!
!  Discussion:
!
!    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
!    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
!    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
!    other abscissas.
!
!    A formal representation is:
!
!      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
!       ( T - T(J) ) / ( T(I) - T(J) )
!
!    This routine accepts a set of INTERP_NUM values at which all the Lagrange
!    polynomials should be evaluated.
!
!    Given data values P_DATA at each of the abscissas, the value of the
!    Lagrange interpolating polynomial at each of the interpolation points
!    is then simple to compute by matrix multiplication:
!
!      P_INTERP(1:INTERP_NUM) =
!        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
!
!    or, in the case where P is multidimensional:
!
!      P_INTERP(1:M,1:INTERP_NUM) =
!        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!    DATA_NUM must be at least 1.
!
!    Input, real ( kind = 8 ) T_DATA(DATA_NUM), the data points.
!
!    Input, integer ( kind = 4 ) INTERP_NUM, the number of
!    interpolation points.
!
!    Input, real ( kind = 8 ) T_INTERP(INTERP_NUM), the
!    interpolation points.
!
!    Output, real ( kind = 8 ) L_INTERP(DATA_NUM,INTERP_NUM), the values
!    of the Lagrange polynomials at the interpolation points.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lagrangeval( xx_value,l_interp,nxyvalue,xold_vec,nold_vec )
    use ModDataType
    implicit none

    integer  (kind=IB), intent(in) :: nold_vec,nxyvalue
    integer  (kind=IB) :: i,j
    
    real     (kind=RP), dimension(nold_vec), intent(in) :: xold_vec
    
    real     (kind=RP), dimension(nxyvalue), intent(in) :: xx_value
    real     (kind=RP), dimension(nold_vec,nxyvalue), intent(out) :: l_interp

    !f2py real     (kind=RP), intent(in)  :: xold_vec
    !f2py real     (kind=RP), intent(in)  :: xx_value
    !f2py real     (kind=RP), intent(out) :: l_interp
    !f2py                     intent(hiden), depend(xold_vec) :: nold_vec=shape(xold_vec,0)
    !f2py                     intent(hiden), depend(xx_value) :: nxyvalue=shape(xx_value,0)
    !f2py                     intent(hiden), depend(l_interp) :: nold_vec=shape(yy_value,0),nxyvalue=shape(yy_value,1)

! *** Evaluate the Lagrange polynomials *************************************
    l_interp(1:nold_vec,1:nxyvalue) = 1.0_RP
    
    do i=1,nold_vec
       do j=1,nold_vec
          
          if ( j /= i ) then
             l_interp(i,1:nxyvalue) = l_interp(i,1:nxyvalue)                &
         * ( xx_value(1:nxyvalue)-xold_vec(j) ) / ( xold_vec(i)-xold_vec(j) )
          end if
          
       end do
    end do
! *** Evaluate the Lagrange polynomials *************************************

    return
END SUBROUTINE lagrangeval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ncc_abscissas ( n, x )

!*****************************************************************************80
!
!! NCC_ABSCISSAS computes the Newton Cotes Closed abscissas.
!
!  Discussion:
!
!    The interval is [ -1, 1 ].
!
!    The abscissas are the equally spaced points between -1 and 1,
!    including the endpoints.
!
!    If N is 1, however, the single abscissas is X = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_ABSCISSAS - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = 0.0D+00
    return
  end if

  do i = 1, n
    x(i) = ( real ( n - i,     kind = 8 ) * ( -1.0D+00 )   &
           + real (     i - 1, kind = 8 ) * ( +1.0D+00 ) ) &
           / real ( n     - 1, kind = 8 )
  end do

  return
end
subroutine ncc_abscissas_ab ( a, b, n, x )

!*****************************************************************************80
!
!! NCC_ABSCISSAS_AB computes the Newton Cotes Closed abscissas for [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCC_ABSCISSAS_AB - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  if ( n == 1 ) then
    x(1) = 0.5D+00 * ( b + a )
    return
  end if

  do i = 1, n
    x(i) = ( real ( n - i,     kind = 8 ) * a   &
           + real (     i - 1, kind = 8 ) * b ) &
           / real ( n     - 1, kind = 8 )
  end do

  return
end
subroutine nco_abscissas ( n, x )

!*****************************************************************************80
!
!! NCO_ABSCISSAS computes the Newton Cotes Open abscissas.
!
!  Discussion:
!
!    The interval is [ -1, 1 ].
!
!    The abscissas are the equally spaced points between -1 and 1,
!    not including the endpoints.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCO_ABSCISSAS - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * ( -1.0D+00 )   &
           + real (     i,     kind = 8 ) * ( +1.0D+00 ) ) &
           / real ( n     + 1, kind = 8 )
  end do

  return
end
subroutine nco_abscissas_ab ( a, b, n, x )

!*****************************************************************************80
!
!! NCO_ABSCISSAS_AB computes the Newton Cotes Open abscissas for [A,B].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, integer ( kind = 4 ) N, the order of the rule.
!
!    Output, real ( kind = 8 ) X(N), the abscissas.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'NCO_ABSCISSAS_AB - Fatal error!'
    write ( *, '(a)' ) '  N < 1.'
    stop 1
  end if

  do i = 1, n
    x(i) = ( real ( n - i + 1, kind = 8 ) * a   &
           + real (     i,     kind = 8 ) * b ) &
           / real ( n     + 1, kind = 8 )
  end do

  return
end
subroutine parameterize_arc_length ( m, data_num, p_data, t_data )

!*****************************************************************************80
!
!! PARAMETERIZE_ARC_LENGTH parameterizes data by pseudo-arclength.
!
!  Discussion:
!
!    A parameterization is required for the interpolation.
!
!    This routine provides a parameterization by computing the
!    pseudo-arclength of the data, that is, the Euclidean distance
!    between successive points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the data values.
!
!    Output, real ( kind = 8 ) T_DATA(DATA_NUM), parameter values
!    assigned to the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) m

  integer ( kind = 4 ) j
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ) t_data(data_num)

  t_data(1) = 0.0D+00
  do j = 2, data_num
    t_data(j) = t_data(j-1) &
      + sqrt ( sum ( ( p_data(1:m,j) - p_data(1:m,j-1) )**2 ) )
  end do

  return
end
subroutine parameterize_index ( m, data_num, p_data, t_data )

!*****************************************************************************80
!
!! PARAMETERIZE_INDEX parameterizes data by its index.
!
!  Discussion:
!
!    A parameterization is required for the interpolation.
!
!    This routine provides a naive parameterization by vector index.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the spatial dimension.
!
!    Input, integer ( kind = 4 ) DATA_NUM, the number of data points.
!
!    Input, real ( kind = 8 ) P_DATA(M,DATA_NUM), the data values.
!
!    Output, real ( kind = 8 ) T_DATA(DATA_NUM), parameter values
!    assigned to the data.
!
  implicit none

  integer ( kind = 4 ) data_num
  integer ( kind = 4 ) m

  integer ( kind = 4 ) j
  real ( kind = 8 ) p_data(m,data_num)
  real ( kind = 8 ) t_data(data_num)

  t_data(1) = 0.0D+00
  do j = 2, data_num
    t_data(j) = real ( j - 1, kind = 8 )
  end do

  return
end
subroutine r8mat_expand_linear2 ( m, n, a, m2, n2, a2 )

!*****************************************************************************80
!
!! R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!    In this version of the routine, the expansion is indicated
!    by specifying the dimensions of the expanded array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), a "small" M by N array.
!
!    Input, integer ( kind = 4 ) M2, N2, the number of rows and columns in A2.
!
!    Output, real ( kind = 8 ) A2(M2,N2), the expanded array, which
!    contains an interpolated version of the data in A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) m2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a2(m2,n2)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  real ( kind = 8 ) r
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) s
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2

  do i = 1, m2

    if ( m2 == 1 ) then
      r = 0.5D+00
    else
      r = real ( i - 1, kind = 8 ) &
        / real ( m2 - 1, kind = 8 )
    end if

    i1 = 1 + int ( r * real ( m - 1, kind = 8 ) )
    i2 = i1 + 1

    if ( m < i2 ) then
      i1 = m - 1
      i2 = m
    end if

    r1 = real ( i1 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    r2 = real ( i2 - 1, kind = 8 ) &
       / real ( m - 1, kind = 8 )

    do j = 1, n2

      if ( n2 == 1 ) then
        s = 0.5D+00
      else
        s = real ( j - 1, kind = 8 ) &
          / real ( n2 - 1, kind = 8 )
      end if

      j1 = 1 + int ( s * real ( n - 1, kind = 8 ) )
      j2 = j1 + 1

      if ( n < j2 ) then
        j1 = n - 1
        j2 = n
      end if

      s1 = real ( j1 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      s2 = real ( j2 - 1, kind = 8 ) &
         / real ( n - 1, kind = 8 )

      a2(i,j) = &
        ( ( r2 - r ) * ( s2 - s ) * a(i1,j1) &
        + ( r - r1 ) * ( s2 - s ) * a(i2,j1) &
        + ( r2 - r ) * ( s - s1 ) * a(i1,j2) &
        + ( r - r1 ) * ( s - s1 ) * a(i2,j2) ) &
        / ( ( r2 - r1 ) * ( s2 - s1 ) )

    end do

  end do

  return
end
function r8vec_ascends_strictly ( n, x )

!*****************************************************************************80
!
!! R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!    Notice the effect of entry number 6 in the following results:
!
!      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
!      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
!      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
!
!      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
!      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the size of the array.
!
!    Input, real ( kind = 8 ) X(N), the array to be examined.
!
!    Output, logical R8VEC_ASCENDS_STRICTLY, is TRUE if the
!    entries of X strictly ascend.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical r8vec_ascends_strictly
  real ( kind = 8 ) x(n)

  do i = 1, n - 1
    if ( x(i+1) <= x(i) ) then
      r8vec_ascends_strictly = .false.
      return
    end if
  end do

  r8vec_ascends_strictly = .true.

  return
end
subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end

subroutine r8vec_expand_linear ( n, x, fat, xfat )
!*****************************************************************************80
!
!! R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    This routine copies the old data, and inserts NFAT new values
!    between each pair of old data values.  This would be one way to
!    determine places to evenly sample a curve, given the (unevenly
!    spaced) points at which it was interpolated.
!
!    An R8VEC is an array of double precision real values.
!
!  Example:
!
!    N = 3
!    NFAT = 2
!
!    X(1:N)        = (/ 0.0,           6.0,             7.0 /)
!    XFAT(1:2*3+1) = (/ 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0 /)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of input data values.
!
!    Input, real ( kind = 8 ) X(N), the original data.
!
!    Input, integer ( kind = 4 ) FAT, the number of data values to interpolate
!    between each pair of original data values.
!
!    Output, real ( kind = 8 ) XFAT((N-1)*(FAT+1)+1), the "fattened" data.
!
  implicit none

  integer ( kind = 4 ) fat
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xfat((n-1)*(fat+1)+1)

  k = 0

  do i = 1, n - 1

    k = k + 1
    xfat(k) = x(i)

    do j = 1, fat
      k = k + 1
      xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(i)     &
                + real (       j,     kind = 8 ) * x(i+1) ) &
                / real ( fat     + 1, kind = 8 )
    end do

  end do

  k = k + 1
  xfat(k) = x(n)

  return
end
subroutine r8vec_expand_linear2 ( n, x, before, fat, after, xfat )

!*****************************************************************************80
!
!! R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
!
!  Discussion:
!
!    This routine starts with a vector of data.
!
!    The intent is to "fatten" the data, that is, to insert more points
!    between successive values of the original data.
!
!    There will also be extra points placed BEFORE the first original
!    value and AFTER that last original value.
!
!    The "fattened" data is equally spaced between the original points.
!
!    The BEFORE data uses the spacing of the first original interval,
!    and the AFTER data uses the spacing of the last original interval.
!
!  Example:
!
!    N = 3
!    BEFORE = 3
!    FAT = 2
!    AFTER = 1
!
!    X    = (/                   0.0,           6.0,             7.0       /)
!    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
!              3 "BEFORE's"      Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 December 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of input data values.
!    N must be at least 2.
!
!    Input, real ( kind = 8 ) X(N), the original data.
!
!    Input, integer ( kind = 4 ) BEFORE, the number of "before" values.
!
!    Input, integer ( kind = 4 ) FAT, the number of data values to interpolate
!    between each pair of original data values.
!
!    Input, integer ( kind = 4 ) AFTER, the number of "after" values.
!
!    Output, real ( kind = 8 ) XFAT(BEFORE+(N-1)*(FAT+1)+1+AFTER), the
!    "fattened" data.
!
  implicit none

  integer ( kind = 4 ) after
  integer ( kind = 4 ) before
  integer ( kind = 4 ) fat
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xfat(before+(n-1)*(fat+1)+1+after)

  k = 0
!
!  Points BEFORE.
!
  do j = 1 - before + fat, fat
    k = k + 1
    xfat(k) = ( real ( fat - j + 1, kind = 8 ) * ( x(1) - ( x(2) - x(1) ) ) &
              + real (       j,     kind = 8 ) *   x(1)          ) &
              / real ( fat     + 1, kind = 8 )
  end do
!
!  Original points and FAT points.
!
  do i = 1, n - 1

    k = k + 1
    xfat(k) = x(i)

    do j = 1, fat
      k = k + 1
      xfat(k) = ( real ( fat - j + 1, kind = 8 ) * x(i)     &
                + real (       j,     kind = 8 ) * x(i+1) ) &
                / real ( fat     + 1, kind = 8 )
    end do

  end do

  k = k + 1
  xfat(k) = x(n)
!
!  Points AFTER.
!
  do j = 1, after
    k = k + 1
    xfat(k) = &
      ( real ( fat - j + 1, kind = 8 ) * x(n)     &
      + real (       j,     kind = 8 ) * ( x(n) + ( x(n) - x(n-1) ) ) ) &
      / real ( fat     + 1, kind = 8 )
  end do

  return
end

FUNCTION r8vec_sorted_nearest ( n, a, value )

!*****************************************************************************80
!
!! R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 April 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Input, real ( kind = 8 ) A(N), a sorted vector.
!
!    Input, real ( kind = 8 ) VALUE, the value whose nearest vector
!    entry is sought.
!
!    Output, integer ( kind = 4 ) R8VEC_SORTED_NEAREST, the index of the nearest
!    entry in the vector.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) r8vec_sorted_nearest
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) lo
  integer ( kind = 4 ) mid
  real ( kind = 8 ) value

  if ( n < 1 ) then
    r8vec_sorted_nearest = -1
    return
  end if

  if ( n == 1 ) then
    r8vec_sorted_nearest = 1
    return
  end if

  if ( a(1) < a(n) ) then

    if ( value < a(1) ) then
      r8vec_sorted_nearest = 1
      return
    else if ( a(n) < value ) then
      r8vec_sorted_nearest = n
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = 1
    hi = n

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return
!
!  A descending sorted vector A.
!
  else

    if ( value < a(n) ) then
      r8vec_sorted_nearest = n
      return
    else if ( a(1) < value ) then
      r8vec_sorted_nearest = 1
      return
    end if
!
!  Seek an interval containing the value.
!
    lo = n
    hi = 1

    do while ( lo < hi - 1 )

      mid = ( lo + hi ) / 2

      if ( value == a(mid) ) then
        r8vec_sorted_nearest = mid
        return
      else if ( value < a(mid) ) then
        hi = mid
      else
        lo = mid
      end if

    end do
!
!  Take the nearest.
!
    if ( abs ( value - a(lo) ) < abs ( value - a(hi) ) ) then
      r8vec_sorted_nearest = lo
    else
      r8vec_sorted_nearest = hi
    end if

    return

  end if

  return
END FUNCTION r8vec_sorted_nearest

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
subroutine timestamp ( )

  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
