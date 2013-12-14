program driver
use base_types, only: dp
use hydro_1d, only: advect
use inits, only: hydro_problem, init
implicit none

integer :: i
type(hydro_problem) :: hp
character*100 :: init_file

! Set up initial conditions
CALL getarg(1, init_file)
call init(hp,init_file)

! Do time update
do i = 1,300
  call advect(hp%x,hp%q,hp%u,hp%dt,hp%flux_lim,hp%bcl,hp%bcr)
end do

! Create new directory and write output
call system('mkdir -p outputs')
open(unit=10, file="outputs/output.txt")
do i=1,size(hp%x)
   write(10,'(2(ES13.4))')  hp%x(i), hp%q(i)
end do
close(10)

end program driver

