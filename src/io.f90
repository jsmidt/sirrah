module io
use base_types, only: dp
implicit none


type hydro_problem

integer :: N 
real(dp), allocatable, dimension(:) :: x, u, q      ! Grid and velocity
integer :: flux_lim, bcl, bcr, dtdump
real(dp) :: dt                      ! Time Step
real(dp) :: tmin, tmax
character*100 :: pname

end type

contains

! Get initial conditions.
subroutine init(this,init_file)
use lib_array, only:  linspace

type(hydro_problem), intent(inout) :: this
character*100, intent(in) :: init_file
real(dp) :: xmin, xmax, dt, tmin, tmax
integer :: N, flux_lim, bcl, bcr, i, dtdump
character*100 :: profile,pname

! Namelist in initalize problem
namelist /input/pname,N,tmin,tmax,dt,dtdump,flux_lim,bcl,bcr,profile

! Read in the namelist from the init file
open(13,file=trim(init_file),status='old')
read(13,nml=input)
close(13)

! Define initial conditions
this%N = N
allocate(this%x(N))
allocate(this%u(N))
allocate(this%q(N))
this%tmin = tmin
this%tmax = tmax
this%dt = dt
this%flux_lim = flux_lim
this%bcl = bcl
this%bcr = bcr
this%dtdump = dtdump
this%pname = pname

! Read in initial profile
open(14,file=trim(profile),status='old')
do i = 1,N
   read(14,*) this%x(i), this%u(i), this%q(i)
end do
close(14)

end subroutine init


! Deallocate arrays
subroutine dalloc_hp(this)
type(hydro_problem), intent(inout) :: this
deallocate (this%x)
deallocate (this%u)
deallocate (this%q)
end subroutine dalloc_hp

subroutine write_hp(this,outfile)
type(hydro_problem), intent(inout) :: this
integer :: i
character*100 :: outfile

where ( abs(this%q) .le. 1d-99 ) this%q = 0.0_dp

! Create new directory and write output
call system('mkdir -p outputs')
open(unit=10, file=outfile)
do i=1,size(this%x)
   write(10,'(2(ES13.4))')  this%x(i), this%q(i)
end do
close(10)

end subroutine write_hp


end module io
