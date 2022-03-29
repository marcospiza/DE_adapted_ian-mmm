!--------------------------------
!..   Fortran .f95 MODULE      ..
!--------------------------------
! MODULE for DE_main.f95
! License: https://github.com/ian-mmm/differential-evolution_f95/blob/master/LICENSE

MODULE moDE !% % % % % % % % %
USE moVarType
IMPLICIT NONE


CONTAINS

!!----------------------------------------------------------
!!..  ftns & subroutines                                  ..
!!----------------------------------------------------------

SUBROUTINE Obj_Ftn(N,BH,val,Vobs,Xin)
! Objective ftn to be MINIMIZED
! EXAMPLE: Griewank function [http://mathworld.wolfram.com/GriewankFunction.html]
! -- requires global vars: dp
    INTEGER,        INTENT(IN)      :: N     ! Number of parameters to be obtimized
    REAL(dp),       INTENT(IN)      :: BH(:) ! Vector with paramter values to be used in the OBJ
    REAL(dp),       INTENT(OUT)     :: val   ! Value returned of the Objective function
    REAL(dp),DIMENSION(:),INTENT(IN):: Vobs  ! Vector with observed values
    REAL(dp),DIMENSION(:),INTENT(IN):: Xin   ! Vector of x-values used to calculate the modelled 

    ! Local variables:
    INTEGER                          :: ii
    REAL(dp)                         :: ymodel               
    REAL(dp)                         :: dummy1,erro
    REAL*8                           :: iRFA
    !--------------------
   
!    if (nxin .ne. nyobs) then
!        print"(A,/,A)","Error in Obj_Ftn subrounine:","size measured and observed values&
!                    &   are different"
!        stop
!    endif 

    PoptVec = BH
    erro = 0.0D0
    DO ii = 1, size(Vobs)
        ! Sub routine to calculate the observed value
        call pb(pvec//"An1",pvec//"Bn1",pvec//"Apm",pvec//"Bpm",&
        &       pvec//"Cpm",pvec//"Dpm",pvec//"coTF",pvec//"alpha",&
        &       pvec//"Kext",pvec//"IAF",Xin(ii),ymodel,dummy1)

        erro = erro +(Vobs(ii) - ymodel  )**2.0D0
    END DO
    val = sqrt(erro)/size(Vobs)

END SUBROUTINE Obj_Ftn
!!======================================================

SUBROUTINE Diff_Evol( nop, NP, LB, UB, T, F_lo, F_hi, Cr, PD, BH_best, F_best,Vobs,Xin)
! Optimizes `functn' using a differential evolution optimization technique
! -- requires global vars: dp
! -- calls ftns/subroutines: Obj_Ftn
    INTEGER,         INTENT(IN)     :: nop
    INTEGER,         INTENT(IN)     :: NP
    REAL(dp),        INTENT(IN)     :: LB(:)
    REAL(dp),        INTENT(IN)     :: UB(:)
    INTEGER,         INTENT(IN)     :: T, PD
    REAL(dp),        INTENT(IN)     :: F_lo, F_hi, Cr
    REAL(dp),        INTENT(OUT)    :: BH_best(:)
    REAL(dp),        INTENT(OUT)    :: F_best
    REAL(dp),DIMENSION(:),INTENT(IN):: Vobs  ! Vector with observed values
    REAL(dp),DIMENSION(:),INTENT(IN):: Xin   ! Vector of x-values used to calculate the modelled 
                                             ! observed value

    ! Local variables:
    INTEGER                                     :: nn, kk, kc, tt,ii
    INTEGER                                     :: IND
    INTEGER,  DIMENSION(1)                      :: f1_pos
    INTEGER,  DIMENSION(3)                      :: spouse
    REAL(dp), DIMENSION( nop, NP )              :: THETA
    REAL(dp), DIMENSION( NP )                   :: F_theta
    REAL(dp), DIMENSION( nop )                  :: range, Z0, theta_new, theta_prime
    REAL(dp), DIMENSION(3)                      :: Z2
    REAL(dp)                                    :: fval, f1, ftri, Z1, RNP, Fdither, &
                                                    RNoP, f1_turn
    character(len=5)                            :: restr
    !--------------------
!   ARGUMENTS:-
!       nop     = (Intent IN) number of total parameters
!       NP      = (Intent IN) number of points in parameter grid
!       LB,UB   = (Intent IN) vector of [lower,upper] bounds for prarameters (beta's)
!       T       = (Intent IN) number of generations
!       PD      = (Intent IN) Indicator for pertubations: 0== off, 1== on
!       SMTH    = (Intent IN) Indicator for using smooth MSE: 0== off, 1== on
!   F_lo, F_hi  = (Intent IN) lower and upper bounds for F dither, "scale factor",
!       Cr      = (Intent IN) crossover value, for binomial where Cr= prob. of spouse/mutant gene
!       eta     = (Intent IN) minimum theshold for ftn value to keep theta vector
!       theta_prime     = "mutant vector" or "spouse vector"
!       theta_new       = "trial vector" or "offspring vector"
!   NOTES:- this diff evol method can be classified as DE/rand/1/bin

RNP     = REAL( NP, dp ) ! convert NP to real
RNoP    = REAL( nop, dp )
f1_turn = 0.0D0
THETA   = 0.0D0

! Create initial candidate solution space:
    ! -- create unchanging variables beforehand outside NP loop
        range = UB - LB
        tt = 0 ! intial grid is generation zero
nn = 1
DO while (nn <= NP) ! -  - @ - - Grid Loop for initial generation - - @ - -
    ! -- parameters uniformly drawn between given bounds
        CALL RANDOM_NUMBER( Z0 )
        THETA( :, nn ) = Z0 * range +  LB
        call subresPm(ResPM,THETA(:, nn))
        if (ResPM) then 
          ! print"(7A10)","An1","Bn1","Apm","Bpm","Cpm","Dpm","Kext"
          ! print"(7F10.5)",THETA(:,nn)

   ! Calculate the function value:
          CALL Obj_Ftn( nop, THETA(:, nn ), fval,Vobs,Xin)
         ! print"(7F10.5)",pvec//"An1",pvec//"Bn1",pvec//"Apm",&
         !               &pvec//"Bpm",pvec//"Cpm",pvec//"Dpm",&
         !               &pvec//"Kext"
          F_theta( nn )   = fval
         !print*,"fval : ",fval, ResPM
          nn = nn + 1
        endif
END DO ! -  - @ - - end grid Loop for initial gen - - @ - -

F_best  = MINVAL( F_theta(:) )

gen_do: DO tt = 1, T !- - generation loop - - - - - - - - - - - - - - - - - - - -
    CALL RANDOM_NUMBER( Z1 )
    Fdither = Z1 * (F_hi - F_lo) + F_lo

    nn = 1
    np_do: DO WHILE( nn <= NP) !~ ~ ~ ~ ~ ~ NP loop ~ ~ ~ ~ ~ ~
        IND = 0
        ! Need to create a SPOUSE for nn:
            ! -- sample from {1,...,NP} without replacement (already removing nn)
            CALL RANDOM_NUMBER( Z2 )
            spouse = FLOOR( Z2 * RNP + 1.0D0)
            DO WHILE (( spouse(1) == nn ))
                CALL RANDOM_NUMBER( Z1 )
                spouse(1) = FLOOR( Z1 * RNP + 1.0D0)
            END DO
            DO WHILE (( spouse(2) == nn ) .AND. ( spouse(2) == spouse(1)))
                CALL RANDOM_NUMBER( Z1 )
                spouse(2) = FLOOR( Z1 * RNP + 1.0D0)
            END DO
            DO WHILE (( spouse(3) == nn ) .AND. ( spouse(3) == spouse(1)) .AND. ( spouse(3) == spouse(2)))
                CALL RANDOM_NUMBER( Z1 )
                spouse(3) = FLOOR( Z1 * RNP + 1.0D0)
            END DO
            ! Standard spouse creation,
                theta_prime(1: nop ) = THETA(1: nop, spouse(1)) + Fdither *( THETA(1: nop, spouse(2)) - THETA(1: nop, spouse(3)) )

        ! Gene transferring:
            CALL RANDOM_NUMBER( Z0 )
            kc = FLOOR( Z0( nop ) * ( nop + 1.0D0) + 1.0D0 ) ! randomly choose first gene to have from spouse/mutatant
            DO kk = 1, nop
                IF ( Z0( kk ) <= Cr .OR. kc == kk ) THEN
                    theta_new( kk ) = theta_prime( kk )
                    ! Now we must test if the new parameter value violates boundary conditions
                    if (theta_new( kk ) > UB ( kk ) .or. theta_new(kk) < LB(kk)) then
                        call random_number( Z1 )
                        theta_new(kk) = LB (kk) + (UB (kk) - LB(kk)) * Z1
                    endif
                ELSE
                    theta_new( kk ) = THETA( kk, nn )
                END IF
            END DO

       ! Test the if the new parameter set satisfies the restriction contions
            call subresPm(ResPM,theta_new(:))
            if (ResPM) then 


        ! Calculate the function value:
               CALL Obj_Ftn( nop, theta_new(:), ftri,Vobs,Xin)
                   ! -- SUBROUTINE Obj_Ftn( N, BH, val)
          !     print*,PoptVec
           ! Picking the Winner:
               IF ( ftri < F_theta( nn ) ) THEN
                   THETA( :, nn ) = theta_new
                   F_theta( nn ) = ftri
                   IND = 1 ! record replacements
               END IF

           ! Recording stats on all vectors,
               f1_turn = f1_turn + REAL(IND,dp)
               nn = nn + 1
          !  else
          !    print*,"Out of restriction"
            endif

    END DO np_do !~ ~ ~ ~ ~ ~ end NP loop ~ ~ ~ ~ ~ ~
    
    ! Calculate averages and grid stats,
        f1              = MINVAL( F_theta ) ! ftn value of best theta for each generation
        f1_pos          = MINLOC( F_theta ) ! location of best vector
        f1_turn         = f1_turn / RNP

 
    ! Live reporting
        PRINT*, tt, f1, f1_pos, f1_turn

!       stop
IF ( f1 < F_best ) THEN
    BH_best(:)  = THETA(:, f1_pos(1))
    F_best  = f1
END IF

! = = = = GRID LOCK EXIT = = = =
IF( tt> 250 ) THEN
    IF (f1_turn <= 0.0000000000000001D0 .AND. PD==0) THEN
        EXIT
    END IF
END IF

    IF (PD ==1) THEN != = = = PERTURBATIONS = = = =
        kk = NINT( REAL(tt,dp) / 10.0D0) * 10
        IF ( tt == kk)  THEN
            IF ( f1_turn < 0.010D0 ) THEN
                PRINT*, "Pertubations! -- -- -- -- -- --"
                DO nn =1, NP !---
                    IF ( nn /= f1_pos(1) ) THEN ! don't throw away the best
                        CALL RANDOM_NUMBER( Z1 )
                        IF ( Z1 > 0.1D0 ) THEN ! chance to keep some bad ones
                            ResPM = .FALSE. 
                            ii = 0
                            do while(.not. ResPM )
                               CALL RANDOM_NUMBER( Z0 )
                               THETA( :, nn ) = Z0 * range + LB
                               ii = ii + 1
                               print*,ii
                               call subresPm(ResPM,THETA(:, nn))
                            enddo

                            CALL Obj_Ftn( nop, THETA(:, nn ), fval,Vobs,Xin)
                            ! -- SUBROUTINE Obj_Ftn( N, BH, val)
                            F_theta( nn ) = fval

                        END IF
                    END IF
                END DO !---------
            END IF
        END IF
    END IF

END DO gen_do !- - - - - - - - end generation loop  - - - - - - - - - - - - - - - - - - - - - - - - -

END SUBROUTINE Diff_Evol
!======================================================

END MODULE moDE !% % % % % % % % %
