subroutine pb(An1,Bn1,Apm,Bpm,Cpm,Dpm,coTF,alpha,Kext,IAF,RFA,fotosynD,Nl)
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

end subroutine pb

