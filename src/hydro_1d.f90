module hydro_1d
use base_types, only: dp
implicit none

contains

! Main flux-limited advection routine
! Inspired by Dullemond's lectures:
! http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/index.shtml
! Equations from these lectures in comments
subroutine advect(x,q,u,dt,flux_lim)
real(dp), dimension(:), intent(in) :: x, u      ! Grid and velocity
integer, intent(in) :: flux_lim
real(dp), dimension(:), intent(inout) :: q      ! Q value
real(dp), dimension(size(q)) :: q_new           ! Temporary new grid
real(dp), dimension(size(u)+2) :: u_half        ! Velocity at interfaces
real(dp), intent(in) :: dt                      ! Time Step
real(dp) :: xp, xm, dx, fluxp, fluxm, r, theta, phir,up, um
integer :: i

! Initialize some things
dx = x(2) - x(1)
q_new = q

! Do the general advection update across the grid.
do i = 3,size(x)-2

    ! Get velcoity at the upper (up) and lower (um) face.
    up = 0.5*(u(i+1)+u(i))
    um = 0.5*(u(i)+u(i-1))

    ! Get x at upper and lower face and form dx. See eq. (4.10)
    xp = 0.5*(x(i+1)+x(i))
    xm = 0.5*(x(i)+x(i-1))
    dx = xp-xm

    ! First find flux at lower interface. (fluxm)
    if (um .ge. 0) then
        theta = 1.0                             ! eq. (4.34)
        r = (q(i-1)-q(i-2))/(q(i)-q(i-1))       ! eq. (4.37)
    else 
        theta = -1.0                            ! eq. (4.34)
        r = (q(i+1)-q(i))/(q(i)-q(i-1))         ! eq. (4.37)
    end if
    phir = get_phir(r,flux_lim)                 ! eq. (4.39-4.41)

    ! Lower flux as seen in eq. (4.38)
    fluxm = 0.5*um*((1+theta)*q(i-1)+(1-theta)*q(i)) &
          + 0.5*abs(um)*(1.0-abs(um*dt/dx))*phir*(q(i)-q(i-1))

    ! First find flux at upper interface. (fluxp)
    if (up .ge. 0) then
        theta = 1.0                             ! eq. (4.34)
        r = (q(i)-q(i-1))/(q(i+1)-q(i))         ! eq. (4.37)
    else 
        theta = -1.0                            ! eq. (4.34)
        r = (q(i+2)-q(i+1))/(q(i+1)-q(i))       ! eq. (4.37)
    end if
    phir = get_phir(r,flux_lim)                 ! eq. (4.39-4.41)

    ! Upper flux as inferred from eq. (4.38)
    fluxp = 0.5*up*((1+theta)*q(i)+(1-theta)*q(i+1)) &
          + 0.5*abs(up)*(1.0-abs(up*dt/dx))*phir*(q(i+1)-q(i))

    ! The general advection update. eq. (4.10)
    q_new(i) = q(i) + dt/dx*(fluxm-fluxp)
end do

! Return the answer.
q = q_new

end subroutine advect


! Choice of a flux limiter
! These flux limiters come from Dullemond's lectures.
! http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/Chapter_4.pdf
! Taken from eq. (4.39-4.41)
real(dp) function get_phir(r,flux_lim)
real(dp), intent(in) :: r
integer :: flux_lim

if (flux_lim .eq. 0) then       ! donor-cell
    get_phir = 0.0
else if (flux_lim .eq. 1) then  ! Lax-Wendroff
    get_phir = 1.0
else if (flux_lim .eq. 2) then  ! Beam-Warming
    get_phir = r
else if (flux_lim .eq. 3) then  ! Fromm
    get_phir = 0.5*(1.0+r)
else if (flux_lim .eq. 4) then  ! minmod
    get_phir = minmod(1.0_dp,r)
else if (flux_lim .eq. 5) then  ! superbee 
    get_phir = max(0.0,min(1.0,2.0*r),min(2.0,r))
else if (flux_lim .eq. 6) then  ! MC
    get_phir = max(0.0,min((1.0+r)/2.0,2.0,2.0*r))
else if (flux_lim .eq. 7) then  ! van Leer
    get_phir = (r+abs(r))/(1.0+abs(r))
else 
    get_phir = 0.0
end if
end function get_phir


! Minmod function. Needed for some flux limiters.
! Equation (4.29)
real(dp) function minmod(a,b)
real(dp), intent(in) :: a, b

if ((abs(a) .lt. abs(b)) .and. (a*b .gt. 0)) then
  minmod = a
else if ((abs(a) .gt. abs(b)) .and. (a*b .gt. 0)) then
  minmod = b
else
  minmod = 0.0_dp
end if
end function minmod


! End of module
end module hydro_1d
