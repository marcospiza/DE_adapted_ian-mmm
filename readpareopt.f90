subroutine readparEopt()
use moVarType
!-------------------------**---------------------------------------------------------------------------
!PURPOSE:
! - Reads the values of each parameter storing into pvec type variable, wich holds name and value
! - Reads the names of the parameters that will need optimization into PoptName char vector
! - Returns the number of parameters that will need optimization into the global integer npOpt
! - Allocates the PoptVec: vector holding the values to be optimized
! - Making pointers between each pvec%val corresponding to the parameter name to be optimized and PoptVec 
!    So that changes in each PoptVec will change the corresponding pvec%val
!-------------------------**---------------------------------------------------------------------------
implicit none

integer:: i ! Local counter
integer:: pos1,pos2 ! Get position to make pointers
integer:: stat


! Open the file form tuttil subroutine
call rdinit(10,20,paramfile)

call readrdval(pvec,"An1",ip)
call readrdval(pvec,"Bn1",ip) 
call readrdval(pvec,"Apm",ip) 
call readrdval(pvec,"Bpm",ip) 
call readrdval(pvec,"Cpm",ip) 
call readrdval(pvec,"Dpm",ip) 
call readrdval(pvec,"CoTF",ip) 
call readrdval(pvec,"alpha",ip) 
call readrdval(pvec,"kext",ip) 
call readrdval(pvec,"IAF",ip) 

! Rreading the parameter names to be optmized
allocate(dummyOpt(maxpar))

call rdacha("PoptName",dummyOpt,maxpar,npOpt)

allocate(PoptName(npOpt))
   PoptName = dummyOpt(1:npOpt)
deallocate(dummyOpt)

! Closing the files from tuttil
close(10) ; close(20)
call RDDTMP(100)


allocate(PoptVec(npOpt))
allocate(BetaLB(npOpt))
allocate(BetaUB(npOpt))
allocate(BetaHAT(npOpt))

! Reading the paroptfile containg set up values for each parameter used in optmizations
call rdinit(10,20,paroptfile)

call readpsetup(pdefV,"An1",isetup)
call readpsetup(pdefV,"Bn1",isetup) 
call readpsetup(pdefV,"Apm",isetup) 
call readpsetup(pdefV,"Bpm",isetup) 
call readpsetup(pdefV,"Cpm",isetup) 
call readpsetup(pdefV,"Dpm",isetup) 
call readpsetup(pdefV,"CoTF",isetup) 
call readpsetup(pdefV,"alpha",isetup) 
call readpsetup(pdefV,"kext",isetup) 
close(10) ; close(20)
call RDDTMP(100)

open(unit=1234, iostat=stat, file="fort.20", status='old')
if (stat == 0) close(1234, status='delete')


do i= 1,npOpt
  ! Making pointers of pvec%val with PoptVec 
  ! Now changing Paramters to te optimized in PoptVec are linked to their respective pvec
   pos1 = posch(PoptName(i),pvec%name) 
   pvec(pos1)%val => PoptVec(i)


  !filling the vectors of lower and upper boundaries of parameter values. Positions os these parameter values
  !are the same as the PoptVec because we are selecting from PoptName sequence
  ! Also setting the initial values for the optmized parameters
  pos2           = posch(PoptName(i),pdefV%name)
  BetaLB(i)      = pdefV(pos2)%vmin
  betaUB(i)      = pdefV(pos2)%vmax
  pvec(pos1)%val = pdefV(pos2)%vi
enddo

end subroutine readparEopt


