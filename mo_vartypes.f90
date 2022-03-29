module moVarType
implicit none

! Precision Parameters: - - - - - - - - - - - - - - - - - - - - - - - - - -
INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 307)
INTEGER, PARAMETER               :: sp = SELECTED_REAL_KIND(6, 37)
                                                            !(P,R) where P=precision and R=decimal exponent range
integer,parameter                :: maxpar = 100            ! Max number of parameters
character(len=14),parameter      :: paroptfile = "paroptsetup.in"
character(len=14),parameter      :: paramfile  = "parameters.dat"


logical                          :: ResPM
integer                          :: ip     = 0               ! counter of parameters, in the parameter file
integer                          :: nobs   = 0               ! counter of parameters, in the input with observed parameters 
integer                          :: isetup = 0               ! counter of parameter names in paroptfile
integer                          :: npOpt                    ! number of paramters to be optmized
character(len=10),allocatable,&                              !Paramter names to be optimized
                & dimension(:)   ::dummyOpt, PoptName        !the dummy is first used
real*8,allocatable,dimension(:),&
                      &    target:: PoptVec                  ! Vector to be optimized. This vector is target to
real*8,allocatable,dimension(:)  :: RFA                      ! PAR for input data
real*8,allocatable,dimension(:)  :: YobsPb                   ! Observed Pb

! ----Variables used in Differential Algorithm method -----------
integer, parameter, dimension(36) :: &
    seedA =(/&
    11981, 31601, 90971, 854099, 295, 4593, &
    369, 12437, 75, 6213, 532, 1989,& 
    11981, 31601, 90971, 854099, 295, 4593, &
    369, 12437, 75, 6213, 532, 1989,& 
    11981, 31601, 90971, 854099, 295, 4593, &
    369, 12437, 75, 6213, 532, 1989/)
INTEGER                           :: T,NP
real(dp)                          :: Fl, Fh, Cr, fval
REAL(dp)                          :: start, finish
real(dp),allocatable,dimension(:) :: BetaUB, BetaLB, BetaHAT


! ******************!!!!**********************************************
!--------------------  Defintion of types -----------------------
type pini
  character(len=10):: name
  real*8:: vi,vmin,vmax                        ! Default values for parameters used in optimization 
end type                                         !vi: initial value; vmin and vmax: minimum and maximum values

type par
  real*8,pointer:: val
  character(len=10):: name
end type

! Type defined variables
type(par),dimension(maxpar):: pvec               ! par type holding name and value for each parameter
type(par),dimension(maxpar):: pobs               ! par type with observed values parameters used for reating yobs
type(pini),dimension(maxpar):: pdefV             ! vector of pini type holding default values from each parameters that
                                                 ! will be used in optmizations

! Creating interfaces assignment for the type par objectcs
interface assignment(=)
 module procedure get_val,val_from_vartype,get_valch,ch_from_vartype
end interface
!-------------------------------------------------------------

! Creating an interface to get the value from the par type vector given the 
! name of the parameter 
interface operator(//)
  module procedure val_from_namepar
end interface
!-------------------------xxx------------------------------------------




contains

! ------- Subroutines for vartype real values assigment -------------------------
subroutine get_val(parType,val)
  ! It sets value to the PartType real value by assigning a value to the detived type
  implicit none 
  type(par),intent(out):: parType
  real*8,intent(in) :: val
  
  if(.not. associated(parType%val)) allocate(parType%val) 
  parType%val = val
end subroutine get_val

subroutine val_from_vartype(val,parType)
! It sets value to a real variable from the derived type
implicit none 
type(par),intent(in):: parType
real*8,intent(out) :: val

val = parType%val
end subroutine val_from_vartype
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

! **** Subroutines for vartype character values assigment ***************8
subroutine get_valch(parType,valch)
! It sets value to the PartType character value by assigning a value to the detived type
implicit none 
type(par),intent(out):: parType
character(len=*),intent(in) :: valch

parType%name = valch
end subroutine get_valch

subroutine ch_from_vartype(ch,parType)
! It sets value to a real variable from the derived type
implicit none 
type(par),intent(in):: parType
character(len=*),intent(out) :: ch

ch = parType%name
end subroutine ch_from_vartype

!*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

!---------------------------------------------------------------------------------------
! Function to be used as new operator to get a value from the var par type by informing
! parameter name 
!---------------------------------------------------------------------------------------

function val_from_namepar(parType,parname)
implicit none
real*8 val_from_namepar
character(len=*),intent(in):: parname
type(par),dimension(:),intent(in):: parType
integer:: pos ! local used to get the position of the parname in the vec

pos = posch(parname,parType%name)
val_from_namepar = parType(pos)%val

end function val_from_namepar

! -------------------------------------------------------------
! Auxiliar functions used in val_from_namepar function
!-------------------------------------------------------------

function posch(val,Vec)
  integer:: posch
  character(*),intent(in):: val
  character(*),dimension(:),intent(in):: vec
  integer:: i
  character(len=len(val)):: varloc

  varloc = to_upper(val)

  i = 0
  do 
   i = i + 1
   if (varloc == to_upper(vec(i)) ) then
      posch = i
      exit
   else
      if (i > size(vec)) then
         print*,"Value not found. You must check ...."
         stop
      endif
   endif
  enddo
end function posch

Pure Function to_upper (str) Result (string)
!   ==============================
!   Changes a string to upper case
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

End Function to_upper

!------------------- ****---------------------------------------
!     End of functions to get out the value form parameter name
!------------------- ****---------------------------------------


! Subrotine to read parameter from the file using ttutil subroutine and save 
! in the parType vector variable
subroutine readrdval(parType,parname,ip)
implicit none
type(par),dimension(:)     :: parType
character(len=*),intent(in):: parname
integer,intent(out)        :: ip    ! Counter of parameters
! ----- Local -------
real*8                     :: parval

ip = ip + 1

parType(ip) = parname
call rdsdou(parname,parval)
parType(ip) = parval

end subroutine readrdval

!************----**************************************************

subroutine readpsetup(parType,parname,isetup)
!------------------------------------------------------------------------------------------------------
! Subrotine to read parameter  from the file containg initial, minimum and maximum value of each 
! parameter used for optimization. The file is a constant-named as paroptsetup.in, defined as parameter 
! by paroptfile
!------------------------------------------------------------------------------------------------------
  implicit none
  type(pini),dimension(:),intent(out):: parType
  character(len=*),intent(in):: parname
  integer,intent(out):: isetup
  
  integer:: ival
  real*8,dimension(3):: dummy
  
  isetup = isetup + 1
  parType(isetup)%name = parname

  call rdfdou(parname,dummy,3,3)
  parType(isetup)%vi   = dummy(1)
  parType(isetup)%vmin = dummy(2)
  parType(isetup)%vmax = dummy(3)
  
end subroutine readpsetup

!************----**************************************************

subroutine readObsY(fname,xin,yobs,ParObs,NpObs)
implicit none
character(len=*)               :: fname
real*8,allocatable,dimension(:):: xin,yobs
type(par),dimension(maxpar)    :: ParObs
integer                        :: NpObs

! Local variables
real*8,dimension(1000):: dummy
integer               :: n1,i


! Open the file form tuttil subroutine
call rdinit(10,20,trim(fname))

! Reading the input and y-observed values
call RDADOU("RFA",dummy,1000,n1)
allocate(xin(n1))
xin = dummy(1:n1)
call RDADOU("PB",dummy,1000,n1)
allocate(yobs(n1))
yobs = dummy(1:n1)

!Converting unities
xin = xin * 1D6


! ****** Raeding parameter values used for y-observed 
call readrdval(ParObs,"An1",NpObs)
call readrdval(ParObs,"Bn1",NpObs) 
call readrdval(ParObs,"Apm",NpObs) 
call readrdval(ParObs,"Bpm",NpObs) 
call readrdval(ParObs,"Cpm",NpObs) 
call readrdval(ParObs,"Dpm",NpObs) 
call readrdval(ParObs,"CoTF",NpObs) 
call readrdval(ParObs,"alpha",NpObs) 
call readrdval(ParObs,"kext",NpObs) 
call readrdval(ParObs,"IAF",NpObs) 

! Closing the files from tuttil
close(10) ; close(20)
call RDDTMP(100)

end subroutine readObsY

subroutine subresPm(res,P_vecOpt)
! Returns the boolean res, which evaluates a restriction about
! the values of the parameters used for computing PM, eq. 9a of LIU qcane paper
logical,intent(out)            :: res
real*8,dimension(:),intent(in) :: P_vecOpt

! Local variables
real*8 :: Rh,Rl

! Assigning the parameters P_vecOpt to the PoptVec

PoptVec = P_vecOpt

Rh = (pvec//"Dpm") * ((pvec//"Bn1") * (pvec//"IAF")) / &
  & ((pvec//"Cpm") * log(1.0 + (pvec//"Bn1") * (pvec//"IAF")))

Rl = (pvec//"Bpm") * ((pvec//"Bn1") * (pvec//"IAF")) / &
  & ((pvec//"Apm") * log(1.0 + (pvec//"Bn1") * (pvec//"IAF")))

res =   Rl < pvec//"An1" .and. pvec//"An1" < Rh

end subroutine subresPm

end module moVarType
