module inits
use base_types, only: dp
implicit none


type hydro_problem

integer :: N 
real(dp), allocatable, dimension(:) :: x, u, q      ! Grid and velocity
integer :: flux_lim
real(dp) :: dt                      ! Time Step
real(dp) :: xmin, xmax, dx

end type

contains

! Get initial conditions.
subroutine init(this,init_file)
use lib_array, only:  linspace

type(hydro_problem), intent(inout) :: this
character*100, intent(in) :: init_file
real(dp) :: xmin, xmax, dt
integer :: N, flux_lim

namelist /input/N,xmin,xmax,dt,flux_lim 

open(13,file=trim(init_file),status='old')
read(13,nml=input)
close(13)

! Define initial conditions
this%N = N

allocate(this%x(this%N))
allocate(this%u(this%N))
allocate(this%q(this%N))

this%xmin = xmin
this%xmax = xmax
call linspace(this%xmin,this%xmax,this%x) ! Set up x-grid
this%q = 0
this%q(1) = 1.0
this%q(2) = 1.0
where ( this%x .le. 30 ) this%q = 1.0_dp
this%u = 1.0
this%dx = this%x(2)-this%x(1)
!dt = 0.1*dx/u(2)
this%dt = dt
this%flux_lim = flux_lim

end subroutine init


end module inits
