! Subrotina para o calculo da fotosintese da cana de acucar
! Programador: Jonathan
! Modified: Marcos Alex at 04-07-20
program photosynt
implicit none
real*8::An1
real*8::Bn1
real*8::IAF
real*8:: RFA(200)
real*8::Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext
real*8:: fotosynD,Nl

integer:: i,nRFA
character(len=10):: obs

! Write output file

open(1,file="simulated.dat",status="replace")
write(1,"(10A)")repeat("-",80)
write(1,*)"     Simulated photosynthesis g CO2/m2 from different PAR"
write(1,"(10A)")repeat("-",80)

open(2,file="observation_pest.dat",status="replace")
write(2,*)"Output for pest pst observation data"

! Read paramter values
call readdata("rfa.dat","parameters.dat",RFA,nRFA,An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF)



do i = 1,nRFA
  call sub_photosynt(An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF,RFA(i)*1e6,fotosynD,Nl)
  
  write(obs,"(I5)")i
  obs = "ar"//trim(adjustl(obs))

  write(1,"(F10.7,2x,F14.7)")RFA(i),fotosynD
  write(2,"(A8,F10.7,f5.1,A8)")obs,fotosynD,1.0,"obs1"
enddo


contains

subroutine sub_photosynt(An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF,RFA,fotosynD,Nl)
implicit none
real*8,intent(in)::An1
real*8,intent(in)::Bn1
real*8,intent(in)::IAF,RFA
real*8,intent(in)::Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext
real*8,intent(out):: fotosynD,Nl

!Local variables
real*8 :: ADco2    !Assimilação Diária de CO2 Bruto por Toda Copa (g CO2 m-2 d-1)
real*8 :: trMM     !Taxa Relativa de MM de CH2O para CO2 (30/40 = 0,682)
real*8 :: Tmn,alphainf,binf
real*8 :: dayleng
real*8 :: Pm,aIAF


dayleng = 12 * 3600 ! converting to seconds
trMM = 0.682    

!redefining alpha
alphainf = 8.28 ! *C
binf = 2.82 ! *C
!Tmn = max(12.0d0,Tmea)
Tmn = 21.0

!alpha = 1.608e-5 * (Tmn -alphainf) / (binf + Tmn)

Nl = (An1 * log(1 + Bn1 * IAF)) / (Bn1 * IAF)

! Correcting with the normalized cilic function 
!Nl = Nl * seascfunc(doy,year(daymeteo),lat,Vsn,Vwn,d1gn) 

! Relação empírica entre Pm x Nl por Alisson et al. (1997) (g CO2 m-2 s-1)
Pm = (Apm * Nl - Bpm) * (1 - exp(Cpm * Nl - Dpm))

! Cálculo do Índice Acumulado de Área foliar
!aIAF = IAF(t) * ((1 - coTF) * Pm * dayleng / (2 * alpha * Kext * RFA(t)))  ! Error in the equation in the paper
aIAF = (1 - coTF) * Pm * dayleng / (2 * alpha * Kext * RFA)

! Cálculo de Assimilação Diária de CO2 Bruto por Toda Copa (g CO2 m-2 d-1)
ADco2  =  (Pm  *  dayleng / Kext) * (Kext * IAF + log((1 + aIAF)/(1 + aIAF * exp(Kext * IAF ))) &
           & + aIAF * log((1 + aIAF) / aIAF) +  aIAF * exp(Kext * IAF)  * &
          &  log(aIAF * exp(Kext * IAF)/(1 + aIAF * exp(Kext * IAF))))

! Calculo de Produçao Diaria de Carboidrato (g CH2O m-2 d-1)
fotosynD  = trMM * ADco2

end subroutine sub_photosynt


subroutine readdata(fname1,fname2,RFA,nRFA,An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF)
implicit none
character(len=*),intent(in):: fname1,fname2
real*8,intent(out):: An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF
real*8,dimension(:),intent(out):: RFA
integer,intent(out):: nRFA

call rdinit(10,20,trim(fname1))
!call rdador("RFA", 0.0d0, 40.0d0, RFA, 200,nRFA)
call rdadou("RFA",RFA,200,nRFA)
close(10)


call rdinit(30,40,trim(fname2))
call rdsdou("An1",An1)
call rdsdou("Bn1",Bn1)
call rdsdou("Apm",Apm)
call rdsdou("Bpm",Bpm)
call rdsdou("Cpm",Cpm)
call rdsdou("Dpm",Dpm)
call rdsdou("coTF",coTF)
call rdsdou("alpha",alpha)
call rdsdou("Kext",Kext)
call rdsdou("IAF",IAF)
close(30)
end subroutine readdata

end
