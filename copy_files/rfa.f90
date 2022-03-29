program rfa_program
implicit none
integer i
real*8 r,RFA

open(1,file = "rfa.dat",status="replace")
write(1,"(2A)")"*",repeat("-",20)
write(1,*) "            RFA =  "
write(1,"(2A)")"*--------- MJ/m2 d",repeat("-",7)

open(2,file="observation_pest.dat",status="replace")

do i = 1,100
  call random_seed()
  call random_number(r)
  r = 10*( r - 0.5)
  RFA =  (10d0 + r)
  write(1,"(F10.7)")RFA
enddo
close(1)

call system("unix2dos rfa.dat")
end
