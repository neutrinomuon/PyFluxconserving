!********************************************************
!*   Program to demonstrate the Akima spline fitting    *
!*       of Function SIN(X) in double precision         *
!* ---------------------------------------------------- *
!*   Reference: BASIC Scientific Subroutines, Vol. II   *
!*   By F.R. Ruckdeschel, BYTE/McGRAWW-HILL, 1981 [1].  *
!*                                                      *
!*                F90 Version By J.-P. Moreau, Paris.   *
!*                        (www.jpmoreau.fr)             *
!* ---------------------------------------------------- *
!* SAMPLE RUN:                                          *
!*                                                      *
!* Akima spline fitting of SIN(X):                      *
!*                                                      *
!*   X   SIN(X) HANDBOOK   AKIMA INTERPOLATION   ERROR  *
!* ---------------------------------------------------- *
!* 0.00    0.0000000          0.0000000       0.0000000 *
!* 0.05    0.0499792          0.0500402      -0.0000610 *
!* 0.10    0.0998334          0.0998435      -0.0000101 *
!* 0.15    0.1494381          0.1494310       0.0000072 *
!* 0.20    0.1986693          0.1986459       0.0000235 *
!* 0.25    0.2474040          0.2474157      -0.0000118 *
!* 0.30    0.2955202          0.2955218      -0.0000016 *
!* 0.35    0.3428978          0.3428916       0.0000062 *
!* 0.40    0.3894183          0.3894265      -0.0000081 *
!* 0.45    0.4349655          0.4349655       0.0000000 *
!* 0.50    0.4794255          0.4794204       0.0000051 *
!* 0.55    0.5226872          0.5226894      -0.0000021 *
!* 0.60    0.5646425          0.5646493      -0.0000068 *
!* 0.65    0.6051864          0.6051821       0.0000043 *
!* 0.70    0.6442177          0.6442141       0.0000035 *
!* 0.75    0.6816388          0.6816405      -0.0000017 *
!* ---------------------------------------------------- *
!*                                                      *
!********************************************************
!PROGRAM Akima
!
!parameter(SIZE=14)
!
!real*8  X(0:SIZE), Y(0:SIZE)
!
!real*8  a,b,xx,yy
!integer i,iv,n
!
!
!  iv=14  ! Number of pooints in table
!  print *,' '
!  print *,'Akima spline fitting of SIN(X):'
!  ! Input sine table
!  !-----------------------------------------------------------------
!  !Sine table values from  Handbook of mathematical functions
!  !by M. Abramowitz and I.A. Stegun, NBS, june 1964
!  !----------------------------------------------------------------- 
!  X(1)=0.000; Y(1)=0.00000000; X(2)=0.125; Y(2)=0.12467473;
!  X(3)=0.217; Y(3)=0.21530095; X(4)=0.299; Y(4)=0.29456472;
!  X(5)=0.376; Y(5)=0.36720285; X(6)=0.450; Y(6)=0.43496553;
!  X(7)=0.520; Y(7)=0.49688014; X(8)=0.589; Y(8)=0.55552980;
!  X(9)=0.656; Y(9)=0.60995199; X(10)=0.721; Y(10)=0.66013615;
!  X(11)=0.7853981634; Y(11)=0.7071067812;
!  X(12)=0.849; Y(12)=0.75062005; X(13)=0.911; Y(13)=0.79011709;
!  X(14)=0.972; Y(14)=0.82601466;
!  !-----------------------------------------------------------------
!  ! print header
!  print *,' '
!  print *,'  X   SIN(X) HANDBOOK   AKIMA INTERPOLATION   ERROR '
!  print *,'----------------------------------------------------'
!  ! main loop
!  xx=0.d0
!  do i = 1, 16
!    call Interpol_Akima(iv,n,xx,yy,X,Y)
!    write(*,50) xx,dsin(xx),yy,dsin(xx)-yy
!    xx = xx + 0.05
!  end do
!  !print footer
!  print *,'----------------------------------------------------'
!  print *,' '
!
!50 format(' ',F4.2,'    ',F9.7,'          ',F9.7,'      ',F10.7)
!  stop
!end


!********************************************************
!*          Akima spline fitting subroutine             *
!* ---------------------------------------------------- *
!* The input table is X(i), Y(i), where Y(i) is the     *
!* dependant variable. The interpolation point is xx,   *
!* which is assumed to be in the interval of the table  *
!* with at least one table value to the left, and three *
!* to the right. The interpolated returned value is yy. *
!* n is returned as an error check (n=0 implies error). *
!* It is also assumed that the X(i) are in ascending    *
!* order.                                               *
!********************************************************
Subroutine Interpol_Akima(NX,n,nxx,xx,yy,X,Y)  
  !Labels: 100,200,300
  INTEGER  :: NX,nxx
  integer i,j
  integer, dimension(0:nxx), intent(out) :: n
  real (kind=8), dimension(0:nxx), intent(in)  :: xx
  real (kind=8), dimension(0:nxx), intent(out) :: yy
  real (kind=8) :: X(0:NX), Y(0:NX)
  real (kind=8) :: XM(0:NX+3)
  real (kind=8) :: Z(0:NX)

  real*8 a,b

  n = 0
  
  do j=0,nxx

     !if ( xx(j) >= X(0) ) then
     !   n(j) = 1
     !   yy(j) = 0.0
     !else if ( xx(j) > X(NX) ) then
     !   n(j) = 1
     !   yy(j) = 0.0
     !else
        !special case xx=0
     !if (xx(j).eq.0.0) then
     !   yy(j)=0.d0 !; return
     !else
     !Check to see if interpolation point is correct
     !if (xx(j)<X(1).or.xx(j)>=X(NX-3)) then
     !   n(j)=0 !; return
     !else !end if
        X(0)=2.d0*X(1)-X(2)
        !Calculate Akima coefficients, a and b
        do i=1,NX!-1
           !Shift i to i+2
           XM(i+2)=(Y(i+1)-Y(i))/(X(i+1)-X(i))
        end do
        XM(NX+2)=2.d0*XM(NX+1)-XM(NX)
        XM(NX+3)=2.d0*XM(NX+2)-XM(NX+1)
        XM(2)=2.d0*XM(3)-XM(4)
        XM(1)=2.d0*XM(2)-XM(3)
        do i=1,NX+1
           a=dabs(XM(i+3)-XM(i+2))
           b=dabs(XM(i+1)-XM(i))
           if (a+b.ne.0.d0) goto 100
           Z(i)=(XM(i+2)+XM(i+1))/2.d0
           goto 200
100        Z(i)=(a*XM(i+1)+b*XM(i+2))/(a+b)
200     end do
        !Find relevant table interval
        i=0
300     i=i+1
        if (xx(j)>X(i)) goto 300
        i=i-1
        !Begin interpolation
        b=X(i+1)-X(i)
        a=xx(j)-X(i)
        yy(j)=Y(i)+Z(i)*a+(3.d0*XM(i+2)-2.d0*Z(i)-Z(i+1))*a*a/b
        yy(j)=yy(j)+(Z(i)+Z(i+1)-2.d0*XM(i+2))*a*a*a/(b*b)
     !end if
  end do
  return

end Subroutine Interpol_Akima
! End of file akima.f90
