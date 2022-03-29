program test2

use moVarType
use moDE

implicit none
integer :: i,pos

real*8    :: test,ymodel,dummy1
character(len=10):: namem


! Fazer subrotina para ler os arquivos de entrada e fala de modo allocatle para quem tenha somente
! o tamanho correto de dados 

! Talvez modificar os vetores para deixarem allocatable. Ler com tamanho maximo com dummy, depois alocar

call readparEopt ! Read parameters and optimization file and more

!----- read input file with observed xin and y-observed and observed parameters
call readObsY("input.dat",RFA,YobsPb,pobs,nobs)



print*,pvec//"An1",pobs//"An1"
print*,(pvec//"kext"),(pobs//"kext")
print*,(pvec//"kext") * (pobs//"kext")




!---------- Optimization procedure ---------------------------
! Start the stopwatch
CALL CPU_TIME( start )

! Optimization Control Variables:
! --Storn and Price (1997), J. of Global Opt., found these rules of thumbs: F \in [0.5, 1.0], Cr \in [0.8, 1.0], Np = 10*D

! Number of population vectors:
NP = 10 * npOpt + 10

! Number of generations
T = 10000

! Crossover variable: Cr \in [0,1]
Cr = 0.85D0 ! manually set
!Cr = 0.85D0 ! manually set

! Scale factor determined by dither:
Fl = SQRT(1.0D0 - 0.50D0* Cr )
Fh = 0.950D0 ! L,H order should be fine as long as Cr > 0.2
! *++*++*++*++*++*++*++*++*++*++*++*

CALL Diff_Evol(npOpt, NP, betaLB, betaUB, T, Fl, Fh, Cr, 1, BetaHat,fval,YobsPb,RFA)
!    CALL Diff_Evol( nop, NP, betaLB, betaUB, T, Fl, Fh, Cr, 1, BetaHAT, fval )
        ! -- SUBROUTINE Diff_Evol( nop, NP, LB, UB, T, F_lo, F_hi, Cr, PD, BH_best, F_best)

PRINT*, "------"
PRINT*, fval, BetaHat
PRINT*, "------"

! Calculate time:
CALL CPU_TIME( finish )
PRINT '("Time = ",f15.2," seconds.")', finish - start

PRINT*, 'DONE ! ! ! ! ! ! ! '



! Results of paramters optimized
open(5,file="de-result.out",status="replace")

write(5,"(A)")"Parameter values measured and optmized"
do i=1,size(PoptName)
   namem = PoptName(i)
   write(5,"(A6, 2E12.3)")namem,pvec//namem,pobs//namem
enddo
write(5,"(A,F18.10)")"OBS = ", fval
! Calculating measured and estimared side by side

write(5,"(A)")"*************************************************************"
write(5,"(A)")"CH2O measured and estimated x RFA "

DO i = 1, size(YobsPb)
    ! Sub routine to calculate the observed value
    call pb(pvec//"An1",pvec//"Bn1",pvec//"Apm",pvec//"Bpm",&
    &       pvec//"Cpm",pvec//"Dpm",pvec//"coTF",pvec//"alpha",&
    &       pvec//"Kext",pvec//"IAF",RFA(i),ymodel,dummy1)
    write(5,"(F10.4,2F12.4)") RFA(i)*1E-6,YobsPb(i),ymodel
END DO

close(5)



!print*,tt(pvec//"kext")



contains
function tt(A)
implicit none
real*8 tt
real*8 A

tt = 2*A

end function tt

end program test2
