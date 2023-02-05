! ###########################################################################
!     RESUME : Rebinning of input vector on a new grid and returns the      !
!              output array. This routine uses the derivative of a          !
!              cumulative function to guarantee the conservation of         !
!              density flux while rebinning the spectrum.                   !
!                                                                           !
!     OBS.  :  Central wavelengths must be specified to compute the         !
!              wavelengths at the edges of the pixels. This routine         !
!              assumes a monotonically increasing function for 'X',         !
!              that can be unequally spaced.                                !
!                                                                           !
!              The Cumulative function has several ways to be               !
!              interpolated:                                                !
!                                                                           !
!              Interpolation Schemes:                                       !
!              00) slow_int -> LINinterpol                                  !
!              01) slow_int -> SPLINE3DFor                                  !
!              02) slow_int -> SPLINE1DArr                                  !
!              03) slow_int -> pchipmodule                                  !
!              04) slow_int -> AkimaSpline                                  !
!              05) slow_int -> Interpolado                                  !
!              06) slow_int -> LINdexerpol                                  !
!              07) slow_int -> poly_interp
!                                                                           !
!     Input           arguments = 7                                         !
!     Output          arguments = 2                                         !
!     Optional        arguments = 2                                         !
!     Total number of arguments = 11                                        !
!                                                                           !
!     INPUT  : 01) Orlambda -> New 'X' vector                               !
!              02) Nrlambda -> Number of elements in new 'X' vector         !
!              03) O_lambda -> Old 'X' vector (Abscissa)                    !
!              04) O_fluxes -> Old 'Y' vector (Ordinate)                    !
!              05) N_lambda -> Number of elements in old 'X' vector         !
!              06) per_bins -> [0: Not conserve flux , 1: Conserve flux  ]  !
!              07) slow_int -> [0 -- 6]                                     !
!              08) fill_val  -> Value to fill values outside boundary       !
!              09) verbosity -> Optional variable to print & check          !
!                                                                           !
!     OUTPUT : 01) Orfluxes -> New 'Y' vector                               !
!              02) IsKeepOn -> Flag, if == 0 then there's a problem         !
!                                                                           !
!     PYTHON : Python compatibility using f2py revised. Better usage        !
!              with numpy.                                                  !
!                                                                           !
!     EXTRA ROUTINES : Those mentioned in slow_int interpolation scheme.    !
!                                                                           !
!     LOG: Modified to take into account Y negative values                  !
!          and minor corrections in the recovered Y values                  !
!          Two different corrections: One related to eps values in the      !
!          cumulative function and another one related to the cumulative    !
!          being set erroneously to zero.                                   !
!                                                                           !
!     Written: Jean Michel Gomes © Copyright ®                              !
!     Checked: Sat Dec  8 12:03:44 WET  2012                                !
!              Fri Dec 28 15:34:27 WET  2012                                !
!              Wed Apr 29 18:13:43 WEST 2015                                !
!              Thu Oct 20 21:55:02 WEST 2016                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FluxConSpec( Orlambda,Orfluxes,Nrlambda,O_lambda,O_fluxes,       &
                        N_lambda,per_bins,slow_int,IsKeepOn,fill_val,       &
                        verbosity )

    use ModDataType
    use pchipmodule
    use Fitting_Fun
    
    implicit none
    integer  (kind=IB), intent(in) :: N_lambda,Nrlambda
    integer  (kind=IB)  :: countdel,indexing,IsShowOn,ilastval,Nz,i
    integer  (kind=IB), optional :: verbosity
    integer  (kind=IB), intent(in) :: per_bins,slow_int
    integer  (kind=IB), intent(out) :: IsKeepOn

    integer  (kind=IB) :: incfd, ierr, findloc_min, findloc_max, Is_Index
    integer  (kind=IB) :: ii, jj, k, yydegree
    integer  (kind=IB), dimension(2) :: next
    
    real     (kind=RP), dimension(N_lambda), intent(in) :: O_lambda,O_fluxes
    real     (kind=RP), dimension(Nrlambda), intent(out) :: Orfluxes
    real     (kind=RP), dimension(Nrlambda), intent(in) :: Orlambda
                                                    
    real     (kind=RP), allocatable, dimension(:) :: auxvecxx,auxvecyy,     &
                                                     auxvecww,auxvecff,     &
                                                     auxvecgg,D_lambda,     &
                                                     Drlambda,O_cumul1,     &
                                                     O_cumul2,Derivada
    real     (kind=RP), optional :: fill_val
    real     (kind=RP) :: max_base,min_base,llow_new,lupp_new,ratioeps,     &
                          sumcumul,maxcumul,tolrance,delta_x,e

    character (len=CH) :: W1aux,W2aux!,string_t
    
    !f2py real     (kind=RP), intent(in) :: Orlambda
    !f2py real     (kind=RP), intent(out) :: Orfluxes
    !f2py                     intent(hide), depend(Orlambda) :: Nrlambda=shape(Orlambda,0)
    !f2py                     intent(hide), depend(Orfluxes) :: Nrlambda=shape(Orfluxes,0)

    !f2py real     (kind=RP), intent(in) :: O_lambda, O_fluxes
    !f2py                     intent(hide), depend(O_lambda) :: N_lambda=shape(O_lambda,0)
    !f2py                     intent(hide), depend(O_fluxes) :: N_lambda=shape(O_fluxes,0)

    !f2py real     (kind=RP), optional :: fill_val=0.0
    
    !f2py integer  (kind=IB), intent(out) :: IsKeepOn
    !f2py integer  (kind=IB), intent(in) :: per_bins,slow_int
    !f2py                     intent(in), optional :: verbosity=0
    
! *** Verbosity mode ********************************************************
    if ( present(verbosity) ) then
        IsShowOn = verbosity
    else
        IsShowOn = 0_IB
    end if
! *** Verbosity mode ********************************************************

    IsKeepOn = 1_IB
    
! *** Check arrays **********************************************************
    !N_fluxes = size(O_fluxes)
    !N_lambda = size(O_lambda)
    !Nrlambda = size(Orlambda)
    !Nrfluxes = size(Orfluxes)

    !if ( N_lambda /= N_fluxes ) then
    !    write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
    !    write (*,'(4x,a)')  '[FluxConSpec] N_fluxes !== N_lambda @@'
    !    write (*,'(4x,a,i10,a,i10)')  '[FluxConSpec]',N_fluxes,' !== ',     &
    !                                  N_lambda
    !
    !    !open  (unit=10,file=arq_outs,position='append',status='unknown')
    !    !write (10,'(4x,a)') '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
    !    !write (10,'(4x,a)') '[FluxConSpec] N_fluxes !== N_lambda @@'
    !    !string_t = fdate()
    !    !write (10,'(4x,a,a)') 'FADO - ',trim(adjustl(string_t))
    !    !close (10)
    !
    !    IsKeepOn = 0
    !    return
    !end if
    !
    !if ( Nrlambda /= Nrfluxes ) then
    !    write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
    !    write (*,'(4x,a)')  '[FluxConSpec] Nrfluxes !== Nrlambda @@'
    !
    !    !open  (unit=10,file=arq_outs,position='append',status='unknown')
    !    !write (10,'(4x,a)') '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
    !    !write (10,'(4x,a)') '[FluxConSpec] Nrfluxes !== Nrlambda @@'
    !    !string_t = fdate()
    !    !write (10,'(4x,a,a)') 'FADO - ',trim(adjustl(string_t))
    !    !close (10)
    !
    !    IsKeepOn = 0
    !    return
    !end if

! *** WARNING ***************************************************************
!     RESUME : This part has been removed to let the code more general.     !
!     N_lambda = count( O_lambda > 0.0_RP )                                 !
!     N_fluxes = N_lambda                                                   !
!     Nrlambda = count( Orlambda > 0.0_RP )                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !write (*,*) N_lambda,Nrlambda
    
    if ( N_lambda < 2_IB .OR. Nrlambda < 2_IB ) then
        write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
        write (*,'(4x,a)')  '[FluxConSpec] not sufficient points @@'
    
        IsKeepOn = 0_IB
        return
    end if

    countdel = count( (O_lambda(2:N_lambda)-O_lambda(1:N_lambda-1))<=0.0_RP )
    if ( countdel > 0_IB ) then
        write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
        write (*,'(4x,a)')  '[FluxConSpec] λ-> Is not monotonically'

        IsKeepOn = 0_IB
        return
    end if
! *** Check arrays **********************************************************

! *** Print screen  *********************************************************
    if ( IsShowOn == 1_IB ) then
        llow_new = Orlambda(1)
        lupp_new = Orlambda(Nrlambda)
        write (*,*)
        write (*,'(4x,a)') '[FluxConSpec]'
        write (W1aux,'(i10)') per_bins
        if ( per_bins == 1 ) then
            write (*,'(4x,a,a,a)') '... per_bins: ',trim(adjustl(W1aux)),   &
                                   ' -- [Conserving flux density]'
        else
            write (*,'(4x,a,a,a)') '... per_bins: ',trim(adjustl(W1aux)),   &
                                   ' -- [Not Conserving flux density]'
        end if
        write (W1aux,'(f17.5)') llow_new
        write (W2aux,'(f17.5)') lupp_new
        write (*,'(4x,5(a))') '[Re-sampling] ==> ',trim(adjustl(W1aux)),    &
                              ' <---> ',trim(adjustl(W2aux)),' Å'
        write (*,'(4x,a,i10,a,i10,a)')                                      &
                     '------------------------------------------------------'
        write (*,'(4x,a,i12,a,i12,a)')                                      &
                     '| N_lambda ==>',N_lambda,' | Nrlambda ==>',Nrlambda,'|'
        write (*,'(4x,a,i10,a,i10,a)')                                      &
                     '------------------------------------------------------'
    end if
! *** Print screen  *********************************************************

! *** Compatibility checks **************************************************
    allocate( auxvecxx(N_lambda+1) )
    allocate( O_cumul1(N_lambda+1) )

    allocate( auxvecyy(Nrlambda+1) )
    allocate( auxvecww(Nrlambda+1) )
    allocate( O_cumul2(Nrlambda+1) )

    allocate( D_lambda(N_lambda) )
    allocate( auxvecff(N_lambda) )

    allocate( Drlambda(Nrlambda) )
    allocate( auxvecgg(Nrlambda) )

! *** Min & Max values ***********
    min_base = max( minval(O_lambda),minval(Orlambda) )
    max_base = min( maxval(O_lambda),maxval(Orlambda) )
    !write (*,*) min_base,max_base

! *** Mid-bin values *************
    D_lambda(2:N_lambda) = O_lambda(2:N_lambda) - O_lambda(1:N_lambda-1)
    D_lambda(1)          = D_lambda(2)
    D_lambda(1:N_lambda) = 0.500000000000000_RP * D_lambda(1:N_lambda)
    auxvecxx(1:N_lambda) = O_lambda(1:N_lambda) - D_lambda(1:N_lambda)
    auxvecxx(N_lambda+1) = O_lambda(N_lambda+0) + D_lambda(N_lambda-1)

! *** Test ******************************************************************
    !do i=1,N_lambda
    !   write (*,*) auxvecxx(i),O_lambda(i),D_lambda(i),i
    !end do
! *** Test ******************************************************************

    if ( count(auxvecxx(2:N_lambda)-auxvecxx(1:N_lambda-1)<=0.0_RP )>0_IB ) then
        write (*,'(4x,a)')  '[PROBLEM_FIT] @@@@@@@@@@@@@@@@@@@@@@@@'
        write (*,'(4x,a)')  '[FluxConSpec] λ is not monotonical @@@'

        IsKeepOn = 0_IB
        return
    end if

    Drlambda(2:Nrlambda) = Orlambda(2:Nrlambda) - Orlambda(1:Nrlambda-1)
    Drlambda(1)          = Drlambda(2)
    Drlambda(1:Nrlambda) = 0.5_RP * Drlambda(1:Nrlambda)
    auxvecyy(1:Nrlambda) = Orlambda(1:Nrlambda) - Drlambda(1:Nrlambda)
    auxvecyy(Nrlambda+1) = Orlambda(Nrlambda+0) + Drlambda(Nrlambda-1)
! *** Compatibility checks **************************************************

! *** Get bin widths ********************************************************
    D_lambda(1:N_lambda) = auxvecxx(2:N_lambda+1) - auxvecxx(1:N_lambda)
    Drlambda(1:Nrlambda) = auxvecyy(2:Nrlambda+1) - auxvecyy(1:Nrlambda)

    if ( per_bins == 1_IB ) then
! *** Unit correction ************
        auxvecff(1:N_lambda) = O_fluxes(1:N_lambda) * D_lambda(1:N_lambda)
    else
        auxvecff(1:N_lambda) = O_fluxes(1:N_lambda)
    end if
! *** Get bin widths ********************************************************

! *** Make cumulative function **********************************************
    tolrance = MachinePrecision( 1_IB )
    O_cumul1(1) = 0.0_RP
    maxcumul    = O_lambda(1)
    sumcumul    = abs( sum( auxvecff(1:N_lambda) ) )
    do indexing=1,N_lambda
        O_cumul1(indexing+1) = O_cumul1(indexing) + auxvecff(indexing)
! --- Test to see whether the original cumulative function varies -----------
        ratioeps = abs(O_cumul1(indexing+1)-O_cumul1(indexing)) / sumcumul
        if ( ratioeps >= tolrance .AND. indexing > 1_IB ) then
           maxcumul = O_lambda(indexing-1)
        end if
! --- Test to see whether the original cumulative function varies -----------
        !write (*,*) O_lambda(indexing),O_cumul1(indexing),sumcumul,ratioeps
    end do
    maxcumul = min( maxcumul,maxval(ARRAY=O_lambda,MASK=O_fluxes/=0.0) )
    maxcumul = int( maxcumul )
    !write (*,*) maxcumul
    !stop
! *** Make cumulative function **********************************************

! *** Interpolate cumulative function ***************************************
    select case (slow_int)
    case default
       ilastval = -999
       call LINinterpol( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,ilastval,IsKeepOn,IsShowOn )
    case (0)
       ilastval = -999
       call LINinterpol( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,ilastval,IsKeepOn,IsShowOn )
    case (1)
       call SPLINE3DArr( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,IsKeepOn,IsShowOn )
    case (2)
       e = 1.0e-8_RP
       Is_Index = 0_IB
       call SPLINE1DArr( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,e,IsKeepOn,IsShowOn )
       
    case (3)
       ! Compute 1st order derivative
       allocate( Derivada(N_lambda+1) )

       incfd = 1
       ierr  = 0
       call dpchim( N_lambda+1, auxvecxx, O_cumul1, Derivada, incfd, ierr )

       next  = 0
       O_cumul2 = 0.0_RP

       do indexing=1,Nrlambda+1
          !findloc_min = minloc( ARRAY=auxvecxx-auxvecyy(indexing), DIM=1,   &
          !                       MASK=auxvecxx>=auxvecyy(indexing) )
          !findloc_max = maxloc( ARRAY=auxvecxx-auxvecyy(indexing), DIM=1,   &
          !                       MASK=auxvecxx<=auxvecyy(indexing) )
          !
          !write (*,*)
          !write (*,*) 'UPPER ==>',findloc_min,findloc_max
          findloc_min = 1
          findloc_max = 0
          
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************
          ii = 1
          if ( auxvecyy(indexing) >= auxvecxx(1) .AND. auxvecyy(indexing) <= auxvecxx(N_lambda+1) ) then
             ii = 1
             jj = N_lambda+1
             do while (jj > ii+1)
                k = (ii+jj)/2
                if (auxvecyy(indexing) < auxvecxx(k)) then
                   jj=k
                else
                   ii=k
                end if
             end do
             findloc_min = ii+1
             findloc_max = ii
          else if ( auxvecyy(indexing) < auxvecxx(0) ) then
             findloc_min = 1
             findloc_max = 0
          else if ( auxvecyy(indexing) > auxvecxx(N_lambda) ) then
             findloc_max = N_lambda+1
             findloc_min = 0 
          end if

          !if ( findloc_min < findloc_max .AND. auxvecyy(indexing) < 7.3 ) then
          !   write (*,*) 'LOWER ==>',findloc_min,findloc_max,indexing,auxvecyy(indexing),auxvecxx(N_lambda+1)
          !end if
          !return
! *** Binary search for i, such that x(i) <= u <= x(i+1) ********************


          !write (*,*) findloc_min, findloc_max, indexing!, !, Orlambda(indexing) 
   
          !Small problem with the indexes and setting the cumulative function to max or zero
          if ( (findloc_min > 0) .AND. (findloc_max > 0) .AND.            &
               (findloc_min /= findloc_max) ) then
             call dchfev( auxvecxx(findloc_min), auxvecxx(findloc_max),     &
                          O_cumul1(findloc_min), O_cumul1(findloc_max),     &
                          Derivada(findloc_min), Derivada(findloc_max),1_IB,&
                          auxvecyy(indexing), O_cumul2(indexing), next, ierr)
          else if ( (findloc_min > 0) .AND. (findloc_max > 0) .AND.         &
                    (findloc_min == findloc_max) ) then
             O_cumul2(indexing) = O_cumul1(findloc_min)
          else
             if ( findloc_min < findloc_max ) then
                O_cumul2(indexing) = O_cumul2(indexing-1)
                !O_cumul2(indexing) = maxval(O_cumul1) !0.0_RP
             end if
          end if
          !if (indexing <= Nrlambda ) then
          !   write (*,*) Orlambda(indexing),auxvecyy(indexing),O_cumul2(indexing),Drlambda(indexing)
          !end if
       end do

       deallocate( Derivada )

    case (4)
       delta_x = 0.1
       call AkimaSpline( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,delta_x,IsKeepOn,IsShowOn )

    case (5)
       Is_Index = 0_IB
       call Interpolado( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,IsKeepOn,Is_Index,IsShowOn )

    case (6)
       call LINdexerpol( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,IsKeepOn,IsShowOn )

    case (7)
       yydegree = 11
       call poly_interp( auxvecyy,O_cumul2,Nrlambda+1,auxvecxx,O_cumul1,    &
                         N_lambda+1,yydegree,IsKeepOn,IsShowOn )

    end select
! *** Interpolate cumulative function ***************************************

! *** Test interpolated cumulative function *********************************
!        do indexing=1,Nrlambda+1
!           if (indexing.LE.N_lambda) then
!              write (*,*) auxvecyy(indexing),O_cumul2(indexing),           &
!                          auxvecxx(indexing),O_cumul1(indexing)
!           else
!              write (*,*) auxvecyy(indexing),O_cumul2(indexing),' 0.0 ',   &
!                          ' 0.0 '
!           end if
!        end do
!        stop
! *** Test interpolated cumulative function *********************************

! *** Undo accumulation *****************************************************
    sumcumul = abs(O_cumul1(N_lambda))
    auxvecgg(1:Nrlambda) = O_cumul2(2:Nrlambda+1)-O_cumul2(1:Nrlambda)
    !write (*,*) auxvecgg(1:Nrlambda)
    
    !if ( count( auxvecgg(1:Nrlambda) < 0.0_RP ) > 0 ) then
    !   write (*,*) "Habemus negativo"
    !else
    !   write (*,*) "Non Habemus negativo"
    !end if
    
    do indexing=1,Nrlambda

       ! ratioeps = abs(sumcumul-O_cumul2(indexing))/sumcumul

       !if ( Orlambda(indexing)  < min_base .OR. Orlambda(indexing)  >      &
       !                           max_base .OR. Orlambda(indexing)  >=     &
       !                          maxcumul .OR. ratioeps < tolrance  .OR.  &
       !    Orlambda(indexing) >=                                          &
       !            int(maxval(ARRAY=O_lambda,MASK=O_fluxes/=0.0_RP)) ) then
       !   !auxvecgg(indexing) = 0.0_RP
       !end if

       if ( per_bins == 1_IB ) then
            auxvecgg(indexing) = auxvecgg(indexing) / Drlambda(indexing)
       end if
       if ( Orlambda(indexing) > O_lambda(N_lambda) ) then
          auxvecgg(indexing) = fill_val !0.0_RP
       end if
       if  ( Orlambda(indexing) < O_lambda(1) ) then
          auxvecgg(indexing) = fill_val !0.0_RP
       end if

       !write (*,*) indexing,auxvecgg(indexing),O_cumul2(indexing),        &
       !            Drlambda(indexing),Orlambda(indexing),ratioeps
   end do
   !stop
! *** Undo accumulation *****************************************************

! *** Export result *********************************************************
   Orfluxes(1:Nrlambda) = auxvecgg(1:Nrlambda)
! *** Export result *********************************************************
   
! *** Print screen **********************************************************
    if ( IsShowOn == 2_IB ) then
        Nz = max(N_lambda,Nrlambda)
        do i=1,Nz
!Print screen **********************************************************
            if ( i <= N_lambda .AND. i <= Nrlambda ) then
                write (*,'(4(e15.4))') O_lambda(i),O_fluxes(i),Orlambda(i), &
                                       Orfluxes(i)
            end if
            if ( i  > N_lambda .AND. i <= Nrlambda ) then
                write (*,'(2(a15),2(e15.4))') ' 0.0 ',' 0.0 ',Orlambda(i),  &
                                                              Orfluxes(i)
            end if
            if ( i <= N_lambda .AND. i  > Nrlambda ) then
                write (*,'(2(e15.4),2(a15))') O_lambda(i),O_fluxes(i),      &
                                              ' 0.0 ',' 0.0 '
            end if
        end do
    end if
    if ( IsShowOn == 1_IB ) then
        write (*,'(4x,a)') '[FluxConSpec]'
    end if
! *** Print screen **********************************************************

! *** Deallocate from memory ************************************************
    deallocate( auxvecgg )
    deallocate( auxvecww )
    deallocate( O_cumul1 )
    deallocate( O_cumul2 )
    deallocate( auxvecff )
    deallocate( D_lambda )
    deallocate( Drlambda )
    deallocate( auxvecxx )
    deallocate( auxvecyy )
! *** Deallocate from memory ************************************************

    return

  contains

! ###########################################################################
!     RESUME : This subroutine uses the MINPACK subroutines.                !
!                                                                           !
! *** Function MachinePrecision ******************************************* !
!                                                                           !
!     This function provides real machine parameters when the               !
!     appropriate set of data statements is activated (by removing the      !
!     c from column 1) and all other data statements are rendered           !
!     inactive. Most of the parameter values were obtained from the         !
!     corresponding Bell Laboratories Port Library function.                !
!                                                                           !
!     The function statement is:                                            !
!                                                                           !
!     real (kind=RP) :: function MachinePrecision( i )                      !
!                                                                           !
!     where                                                                 !
!                                                                           !
!     i is an integer input variable set to 1, 2, or 3 which selects        !
!     the desired machine parameter. If the machine has t base b            !
!     digits and its smallest and largest exponents are emin and emax,      !
!     respectively, then these parameters are                               !
!                                                                           !
!     MachinePrecision(1) = b**(1 - t), the machine precision,              !
!                                                                           !
!     MachinePrecision(2) = b**(emin - 1), the smallest magnitude,          !
!                                                                           !
!     MachinePrecision(3) = b**emax*(1 - b**(-t)), the largest magnitude.   !
!                                                                           !
! ------------------------------------------------------------------------- !
!                                                                           !
!     Machine constants for the IBM 360/370 series, the Amdahl 470/V6,      !
!     the ICL 2900, the Itel AS/6, the Xerox Sigma 5/7/9 and the Sel        !
!     systems 85/86.                                                        !
!                                                                           !
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /                     !
!     data minmag(1),minmag(2) / z00100000, z00000000 /                     !
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /                     !
!                                                                           !
!     Machine constants for the Honeywell 600/6000 series.                  !
!                                                                           !
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /             !
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /             !
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /             !
!                                                                           !
!     Machine constants for the CDC 6000/7000 series.                       !
!                                                                           !
!     data mcheps(1) / 15614000000000000000b /                              !
!     data mcheps(2) / 15010000000000000000b /                              !
!                                                                           !
!     data minmag(1) / 00604000000000000000b /                              !
!     data minmag(2) / 00000000000000000000b /                              !
!                                                                           !
!     data maxmag(1) / 37767777777777777777b /                              !
!     data maxmag(2) / 37167777777777777777b /                              !
!                                                                           !
!     Machine constants for the PDP-10 (KA processor).                      !
!                                                                           !
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /             !
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /             !
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /             !
!                                                                           !
!     Machine constants for the PDP-10 (KI processor).                      !
!                                                                           !
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /             !
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /             !
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /             !
!                                                                           !
!     Machine constants for the PDP-11.                                     !
!                                                                           !
!     data mcheps(1),mcheps(2) /   9472,      0 /                           !
!     data mcheps(3),mcheps(4) /      0,      0 /                           !
!                                                                           !
!     data minmag(1),minmag(2) /    128,      0 /                           !
!     data minmag(3),minmag(4) /      0,      0 /                           !
!                                                                           !
!     data maxmag(1),maxmag(2) /  32767,     -1 /                           !
!     data maxmag(3),maxmag(4) /     -1,     -1 /                           !
!                                                                           !
!     Machine constants for the Burroughs 6700/7700 systems.                !
!                                                                           !
!     data mcheps(1) / o1451000000000000 /                                  !
!     data mcheps(2) / o0000000000000000 /                                  !
!                                                                           !
!     data minmag(1) / o1771000000000000 /                                  !
!     data minmag(2) / o7770000000000000 /                                  !
!                                                                           !
!     data maxmag(1) / o0777777777777777 /                                  !
!     data maxmag(2) / o7777777777777777 /                                  !
!                                                                           !
!     Machine constants for the Burroughs 5700 system.                      !
!                                                                           !
!     data mcheps(1) / o1451000000000000 /                                  !
!     data mcheps(2) / o0000000000000000 /                                  !
!                                                                           !
!     data minmag(1) / o1771000000000000 /                                  !
!     data minmag(2) / o0000000000000000 /                                  !
!                                                                           !
!     data maxmag(1) / o0777777777777777 /                                  !
!     data maxmag(2) / o0007777777777777 /                                  !
!                                                                           !
!     Machine constants for the Burroughs 1700 system.                      !
!                                                                           !
!     data mcheps(1) / zcc6800000 /                                         !
!     data mcheps(2) / z000000000 /                                         !
!                                                                           !
!     data minmag(1) / zc00800000 /                                         !
!     data minmag(2) / z000000000 /                                         !
!                                                                           !
!     data maxmag(1) / zdffffffff /                                         !
!     data maxmag(2) / zfffffffff /                                         !
!                                                                           !
!     Machine constants for the Univac 1100 series.                         !
!                                                                           !
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /             !
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /             !
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /             !
!                                                                           !
!     Machine constants for the Data General Eclipse S/200.                 !
!                                                                           !
!     Note - it may be appropriate to include the following card -          !
!     static dmach(3)                                                       !
!                                                                           !
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/                         !
!     data mcheps/32020k,3*0/                                               !
!                                                                           !
!     Machine constants for the Harris 220.                                 !
!                                                                           !
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /                     !
!     data minmag(1),minmag(2) / '20000000, '00000201 /                     !
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /                     !
!                                                                           !
!     Machine constants for the Cray-1.                                     !
!                                                                           !
!     data mcheps(1) / 0376424000000000000000b /                            !
!     data mcheps(2) / 0000000000000000000000b /                            !
!                                                                           !
!     data minmag(1) / 0200034000000000000000b /                            !
!     data minmag(2) / 0000000000000000000000b /                            !
!                                                                           !
!     data maxmag(1) / 0577777777777777777777b /                            !
!     data maxmag(2) / 0000007777777777777776b /                            !
!                                                                           !
!     Machine constants for the Prime 400.                                  !
!                                                                           !
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /               !
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /               !
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /               !
!                                                                           !
!     Machine constants for the VAX-11.                                     !
!                                                                           !
!     data mcheps(1),mcheps(2) /   9472,  0 /                               !
!     data minmag(1),minmag(2) /    128,  0 /                               !
!     data maxmag(1),maxmag(2) / -32769, -1 /                               !
!                                                                           !
!     Machine constants for IEEE machines.                                  !
!                                                                           !
!     Written: Argonne national laboratory. minpack project. march 1980.    !
!              Burton s. Garbow, Kenneth E. Hillstrom, Jorge J.             !
!                                                                           !
!     Written: Readapted to f90 by Jean Michel Gomes                        !
!     Checked: Thu Oct 27 09:46:50 WEST 2016                                !
! ###########################################################################
REAL (kind=RP) FUNCTION MachinePrecision( i )
    use ModDataType
    implicit none
    integer  (kind=IB), intent(in) :: i
    integer  (kind=IB) :: mcheps(4)
    integer  (kind=IB) :: minmag(4)
    integer  (kind=IB) :: maxmag(4)
    real     (kind=RP) :: dmach(3)

    !f2py integer  (kind=IB), intent(in) :: i
    !f2py real     (kind=RP), intent(out) :: MachinePrecision
    
    equivalence ( dmach(1),mcheps(1) )
    equivalence ( dmach(2),minmag(1) )
    equivalence ( dmach(3),maxmag(1) )
    
    data dmach(1) /2.22044604926e-016_RP/
    data dmach(2) /2.22507385852e-308_RP/
    data dmach(3) /1.79769313485e+308_RP/
    
    MachinePrecision = dmach(i)
    return
  END FUNCTION MachinePrecision
! ###########################################################################

END SUBROUTINE FluxConSpec
! ###########################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE author_FluxConSpec( a )
  use ModDataType

  implicit none
  
  character (len=21), intent(out) :: a

  !f2py intent(out) :: a

  a = 'Written by Jean Gomes'
  
END SUBROUTINE author_FluxConSpec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Jean@Porto - Wed Dec  5 11:04:34 WET 2012 +++++++++++++++++++++++++++++++++

! *** Test ******************************************************************
!PROGRAM GeneralTest
!END PROGRAM GeneralTest
! *** Test ******************************************************************

! *** Number : 002                                                          !
!
! 1) FluxConSpec
! 2) author_FluxConSpec

