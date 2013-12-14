module hydro_1d
use base_types, only: dp
implicit none

contains

! Main flux-limited advection routine
! Inspired by Dullemond's lectures:
! http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/index.shtml
! Equations from these lectures in comments
subroutine advect(x,q,u,dt,flux_lim,bcl,bcr)
real(dp), dimension(:), intent(in) :: x, u      ! Grid and velocity
integer, intent(in) :: flux_lim, bcl, bcr
real(dp), dimension(:), intent(inout) :: q      ! Q value
real(dp), dimension(size(u)+4) :: q_new, x_gh, q_gh, u_gh   ! Velocity at interfaces
real(dp), intent(in) :: dt                      ! Time Step
real(dp) :: xp, xm, dx, fluxp, fluxm, r, theta, phir,up, um
integer :: i

! Initialize some things
q_new = 0
x_gh = 0
q_gh = 0
u_gh = 0
q_gh(3:size(x_gh)-2) = q
x_gh(3:size(x_gh)-2) = x
u_gh(3:size(x_gh)-2) = u
q_new = q_gh



! Do the general advection update across the grid.
do i = 3,size(x_gh)-2

    ! Update Boundry Conditions
    call update_boudries(x_gh,u_gh,q_gh,bcl,bcr)

    ! Get velcoity at the upper (up) and lower (um) face.
    up = 0.5*(u_gh(i+1)+u_gh(i))
    um = 0.5*(u_gh(i)+u_gh(i-1))

    ! Get x at upper and lower face and form dx. See eq. (4.10)
    xp = 0.5*(x_gh(i+1)+x_gh(i))
    xm = 0.5*(x_gh(i)+x_gh(i-1))
    dx = xp-xm

    ! First find flux at lower interface. (fluxm)
    if (um .ge. 0) then
        theta = 1.0                                     ! eq. (4.34)
        r = (q_gh(i-1)-q_gh(i-2))/(q_gh(i)-q_gh(i-1)+1d-70)   ! eq. (4.37)
    else 
        theta = -1.0                                    ! eq. (4.34)
        r = (q_gh(i+1)-q_gh(i))/(q_gh(i)-q_gh(i-1)+1d-70)     ! eq. (4.37)
    end if
    phir = get_phir(r,flux_lim)                         ! eq. (4.39-4.41)

    ! Lower flux as seen in eq. (4.38)
    fluxm = 0.5*um*((1+theta)*q_gh(i-1)+(1-theta)*q_gh(i)) &
          + 0.5*abs(um)*(1.0-abs(um*dt/dx))*phir*(q_gh(i)-q_gh(i-1))

    ! First find flux at upper interface. (fluxp)
    if (up .ge. 0) then
        theta = 1.0                                     ! eq. (4.34)
        r = (q_gh(i)-q_gh(i-1))/(q_gh(i+1)-q_gh(i)+1d-70)     ! eq. (4.37)
    else 
        theta = -1.0                                    ! eq. (4.34)
        r = (q_gh(i+2)-q_gh(i+1))/(q_gh(i+1)-q_gh(i)+1d-70)   ! eq. (4.37)
    end if
    phir = get_phir(r,flux_lim)                         ! eq. (4.39-4.41)

    ! Upper flux as inferred from eq. (4.38)
    fluxp = 0.5*up*((1+theta)*q_gh(i)+(1-theta)*q_gh(i+1)) &
          + 0.5*abs(up)*(1.0-abs(up*dt/dx))*phir*(q_gh(i+1)-q_gh(i))

    ! The general advection update. eq. (4.10)
    q_new(i) = q_gh(i) + dt/dx*(fluxm-fluxp)
end do

! Return the answer.
q = q_new(3:size(x_gh)-2)

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


! Boundry condition definitions are taken from Dullemond lectures
! section 5.3:
! http://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2012/Chapter_5.pdf
subroutine update_boudries(x_gh,u_gh,q_gh,bcl,bcr)
real(dp), dimension(:), intent(inout) :: x_gh,u_gh,q_gh
integer, intent(in) :: bcl,bcr
integer :: end0,end1,end2,end3,end4
real(dp) :: dx1, dx2

end0 = size(x_gh)
end1 = size(x_gh)-1
end2 = size(x_gh)-2
end3 = size(x_gh)-3
end4 = size(x_gh)-4

! Left Boundry (see section 5.3)
if (bcl .eq. 1) then                            ! Periodic
   dx1 = x_gh(end3) - x_gh(end3) 
   dx2 = x_gh(end3) - x_gh(end4) 
   x_gh(2) = x_gh(3) - dx1
   x_gh(1) = x_gh(2) - dx2
   u_gh(2) = u_gh(end2)
   u_gh(1) = u_gh(end3)
   q_gh(2) = q_gh(end2)
   q_gh(1) = q_gh(end3)
else if (bcl .eq. 2) then                       ! Reflexive
   dx1 = x_gh(4) - x_gh(3) 
   dx2 = x_gh(5) - x_gh(4) 
   x_gh(2) = x_gh(3) - dx1
   x_gh(1) = x_gh(2) - dx2
   u_gh(2) = -u_gh(3)
   u_gh(1) = -u_gh(4)
   q_gh(2) = q_gh(3)
   q_gh(1) = q_gh(4)
else if (bcl .eq. 3) then                     ! Free outflow/inflow
   dx1 = x_gh(4) - x_gh(3) 
   dx2 = x_gh(5) - x_gh(4) 
   x_gh(2) = x_gh(3) - dx1
   x_gh(1) = x_gh(2) - dx2
   u_gh(2) = u_gh(3)
   u_gh(1) = u_gh(4)
   q_gh(2) = q_gh(3)
   q_gh(1) = q_gh(4)
else if (bcl .eq. 4) then                     ! Free outflow, no inflow
   dx1 = x_gh(4) - x_gh(3) 
   dx2 = x_gh(5) - x_gh(4) 
   x_gh(2) = x_gh(3) - dx1
   x_gh(1) = x_gh(2) - dx2
   u_gh(2) = -abs(u_gh(3))
   u_gh(1) = -abs(u_gh(4))
   q_gh(2) = q_gh(3)
   q_gh(1) = q_gh(4)
end if

! Right Boundry (see section 5.3)
if (bcr .eq. 1) then                            ! Periodic
   dx1 = x_gh(4) - x_gh(3)
   dx2 = x_gh(5) - x_gh(4)
   x_gh(end1) = x_gh(end2) + dx1
   x_gh(end0) = x_gh(end1) + dx2
   u_gh(end1) = u_gh(3)
   u_gh(end0) = u_gh(4)
   q_gh(end1) = q_gh(3)
   q_gh(end0) = q_gh(4)
else if (bcr .eq. 2) then                     ! Reflexive
   dx1 = x_gh(end2) - x_gh(end3)      
   dx2 = x_gh(end3) - x_gh(end4)      
   x_gh(end1) = x_gh(end2) + dx1
   x_gh(end0) = x_gh(end1) + dx2
   u_gh(end1) = -u_gh(end2)
   u_gh(end0) = -u_gh(end3)
   q_gh(end1) = q_gh(end2)
   q_gh(end0) = q_gh(end3)
else if (bcr .eq. 3) then                     ! Free outflow/inflow
   dx1 = x_gh(end2) - x_gh(end3)      
   dx2 = x_gh(end3) - x_gh(end4)      
   x_gh(end1) = x_gh(end2) + dx1
   x_gh(end0) = x_gh(end1) + dx2
   u_gh(end1) = u_gh(end2)
   u_gh(end0) = u_gh(end3)
   q_gh(end1) = q_gh(end2)
   q_gh(end0) = q_gh(end3)
else if (bcr .eq. 4) then                     ! Free outflow, no inflow
   dx1 = x_gh(end2) - x_gh(end3)      
   dx2 = x_gh(end3) - x_gh(end4)      
   x_gh(end1) = x_gh(end2) + dx1
   x_gh(end0) = x_gh(end1) + dx2
   u_gh(end1) = -abs(u_gh(end2))
   u_gh(end0) = -abs(u_gh(end3))
   q_gh(end1) = q_gh(end2)
   q_gh(end0) = q_gh(end3)
end if

end subroutine update_boudries


! End of module
end module hydro_1d
