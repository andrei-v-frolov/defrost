! $Id: defrost.f90,v 1.1 2007/04/30 23:21:24 frolov Exp $
! [compile with: ifc -O3 -ipo -pc80 -r8 -w95 defrost.f90]

! Reheating code doing something...

program main; implicit none

! control parameters
integer, parameter :: n	= 256			! grid size: spatial grid contains n^3 points
integer, parameter :: nt = 32			! grid size: spatial grid contains n^3 points
integer, parameter :: nx = 32			! grid size: spatial grid contains n^3 points
integer, parameter :: p	= n+2			! grid padding: must be
real, parameter :: alpha = 2.0
real, parameter :: dx = 1.0/n, dt = dx/alpha

integer l
integer i, j, k

! smp array samples all fields on a 3D grid for three subsequent time slices
integer, parameter :: phi = 1, psi = 2, fields = 2; real smp(fields,0:p,0:p,0:p,3)

! for larger grids, dynamically allocated array might be required
! real, allocatable :: smp(:,:,:,:,:); allocate(smp(fields,0:p,0:p,0:p,3))

smp = 0.0
smp(1,n/2,n/2,n/2,:) = 1.0

call random_number(smp(:,:,:,:,1))

!smp(1,:,:,:,2) = sqrt(-2.0*log(smp(1,:,:,:,1)))
!smp(2,:,:,:,2) = smp(1,:,:,:,2)*sin(smp(2,:,:,:,1))
!smp(1,:,:,:,2) = smp(1,:,:,:,2)*cos(smp(2,:,:,:,1))

do k = 1,n; do j = 1,n; do i = 1,n
	smp(:,i,j,k,2) = sqrt(-2.0*log(smp(1,i,j,k,1))) * (/ sin(smp(2,i,j,k,1)), cos(smp(2,i,j,k,1)) /)
end do; end do; end do

do l = 1,10,3
	!call step(l,   smp(:,:,:,:,1), smp(:,:,:,:,2), smp(:,:,:,:,3))
	!call step(l+1, smp(:,:,:,:,2), smp(:,:,:,:,3), smp(:,:,:,:,1))
	!call step(l+2, smp(:,:,:,:,3), smp(:,:,:,:,1), smp(:,:,:,:,2))
end do

contains

! scalar field evolution step
subroutine step(l, dn, hr, up)
	real, dimension(fields,0:p,0:p,0:p) :: dn, hr, up; integer i, j, k, l
	
	real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -6.0, cc = 1.0
	!real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 2.0, c0 = -24.0, cc = 6.0
	!real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 8.0, c0 = -56.0, cc = 12.0
	!real, parameter :: c3 = 1.0, c2 = 3.0, c1 = 14.0, c0 = -128.0, cc = 30.0
	
	real, parameter :: a = c0 + 2.0 * alpha**2 * cc, c = alpha**2 * cc
	
	do k = 1,n; do j = 1,n; do i = 1,n
		up(:,i,j,k) = (                                                                                              &
		    c1 * ( hr(:,i-1,j,k) + hr(:,i+1,j,k) + hr(:,i,j-1,k) + hr(:,i,j+1,k) + hr(:,i,j,k-1) + hr(:,i,j,k+1) ) + &
		    c2 * ( hr(:,i-1,j-1,k) + hr(:,i+1,j-1,k) + hr(:,i-1,j+1,k) + hr(:,i+1,j+1,k)                             &
		         + hr(:,i-1,j,k-1) + hr(:,i+1,j,k-1) + hr(:,i-1,j,k+1) + hr(:,i+1,j,k+1)                             &
			 + hr(:,i,j-1,k-1) + hr(:,i,j+1,k-1) + hr(:,i,j-1,k+1) + hr(:,i,j+1,k+1) ) +                         &
		    c3 * ( hr(:,i-1,j-1,k-1) + hr(:,i+1,j-1,k-1) + hr(:,i-1,j+1,k-1) + hr(:,i+1,j+1,k-1)                     &
		         + hr(:,i-1,j-1,k+1) + hr(:,i+1,j-1,k+1) + hr(:,i-1,j+1,k+1) + hr(:,i+1,j+1,k+1) ) +                 &
		    a * hr(:,i,j,k) )/c - dn(:,i,j,k)
		
		up(phi,i,j,k) = up(phi,i,j,k) - (1.0 + hr(psi,i,j,k)**2) * hr(phi,i,j,k) * dt**2
		up(psi,i,j,k) = up(psi,i,j,k) - (2.0 + hr(phi,i,j,k)**2) * hr(psi,i,j,k) * dt**2
	end do; end do; end do
	
	up(:,0,:,:) = up(:,n,:,:); up(:,n+1,:,:) = up(:,1,:,:)
	up(:,:,0,:) = up(:,:,n,:); up(:,:,n+1,:) = up(:,:,1,:)
	up(:,:,:,0) = up(:,:,:,n); up(:,:,:,n+1) = up(:,:,:,1)
	
	! ...
	if (mod(l-1, n/nt) == 0) then; call bov(phi, "pHi", l, (l-1)*dt, hr); end if
end subroutine step

! 
subroutine bov(v, name, frame, t, array)
	integer k, v, frame; real t, array(fields,0:p,0:p,0:p); character(*) :: name; character(256) :: buffer
	
	write (buffer,'(a,a,i6.6,a)') name, '-', frame, '.bov'; open(10, file=buffer)
	write (buffer,'(a,a,i6.6,a)') name, '-', frame, '.dat'; open(11, file=buffer, form="binary")
	
	write (10,'(4g)') "TIME:        ", t*dt
	write (10,'(4g)') "BRICK_ORIGIN:", 0.0, 0.0, 0.0
	write (10,'(4g)') "BRICK_SIZE:  ", n*dx, n*dx, n*dx
	write (10,'(4g)') "VARIABLE:        ", name
	write (10,'(4g)') "DATA_FILE:       ", buffer
	write (10,'(g,3i5)') "DATA_SIZE:     ", nx, nx, nx
	write (10,'(4g)') "DATA_FORMAT:     ", "DOUBLE"
	write (10,'(4g)') "DATA_ENDIAN:     ", "LITTLE"
	write (10,'(4g)') "CENTERING:       ", "zonal"
	close (10)
	
	do k = 1,n,n/nx; write (11) array(v,1:n:n/nx,1:n:n/nx,k); end do
	close (11)
end subroutine bov

end
