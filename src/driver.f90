program driver
use base_types, only: dp
use hydro_1d, only: advect
use io, only: hydro_problem, init, dalloc_hp, write_hp
implicit none

integer :: i, dumpn
real(dp) ::  ctime
type(hydro_problem) :: hp
character*100 :: init_file, outfile

! Set up initial conditions
call getarg(1, init_file)
call init(hp,init_file)

! Do time update
ctime = hp%tmin
dumpn = 0
do while (ctime .le. hp%tmax)
  call advect(hp%x,hp%q,hp%u,hp%dt,hp%flux_lim,hp%bcl,hp%bcr)

  if (mod(dumpn,hp%dtdump) .eq. 0) then
     write(outfile,'(A,I6.6,A)') 'outputs/'//trim(hp%pname)//'-',dumpn, '.txt'
     call write_hp(hp,outfile)
  end if 
  ctime = ctime + hp%dt
  dumpn = dumpn + 1
enddo
write(outfile,'(A,I6.6,A)') 'outputs/'//trim(hp%pname)//'-',dumpn, '.txt'
call write_hp(hp,outfile)


call dalloc_hp(hp)

end program driver

