! $Id: defrost.f90,v 1.2 2007/04/30 23:43:43 frolov Exp $
! [compile with: ifc -O3 -ipo -pc80 -r8 -w95 defrost.f90]

! Reheating code doing something...

program main; implicit none

include "fftw3.f"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! control parameters
integer, parameter :: n = 256                   ! grid size: spatial grid contains n^3 points
integer, parameter :: nn = n/2+1, ns = 1.732051*(n/2) + 1
integer, parameter :: nt = 32                   ! grid size: spatial grid contains n^3 points
integer, parameter :: nx = 32                   ! grid size: spatial grid contains n^3 points
integer, parameter :: p = n+2                   ! grid padding: must be
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

do l = 1,10,3
        !call step(l,   smp(:,:,:,:,1), smp(:,:,:,:,2), smp(:,:,:,:,3))
        !call step(l+1, smp(:,:,:,:,2), smp(:,:,:,:,3), smp(:,:,:,:,1))
        !call step(l+2, smp(:,:,:,:,3), smp(:,:,:,:,1), smp(:,:,:,:,2))
end do

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sample Gaussian random field
subroutine sample(f, sigma)
        real f(n,n,n), sigma; complex Fk(nn,n,n); integer*8 plan
        integer i, j, k; real :: s = 0.0, a(nn), p(nn)
        
        complex, parameter :: w = (0.0, 6.2831853071795864769252867665590)
        
        do k = 1,n; do j = 1,n; do i = 1,n
                f(i,j,k) = (16.0 + (i-nn)**2 + (j-nn)**2 + (k-nn)**2)**(-3.0/4.0)
                s = s + f(i,j,k)**2
        end do; end do; end do
        
        s = sigma/sqrt(s*real(n)**3)
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,f,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        call random_seed
        do k = 1,n; do j = 1,n
                call random_number(a); call random_number(p)
                Fk(:,j,k) = s * sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
        end do; end do
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine sample

! solve Laplace equation $\Delta f = \rho$
subroutine laplace(f, rho)
        real f(n,n,n), rho(n,n,n); complex Fk(nn,n,n); integer*8 plan
        integer i, j, k; real :: ii, jj, kk
        
        real, parameter :: pi = 3.1415926535897932384626433832795
        
        real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -3.0, cc = 0.5
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 1.0, c0 = -6.0, cc = 1.5
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 2.0, c0 = -7.0, cc = 1.5
        !real, parameter :: c3 = 2.0, c2 = 3.0, c1 = 7.0, c0 = -32.0, cc = 7.5
        
        real, parameter :: c = cc * dx**2/real(n)**3
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,rho,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        do k = 1, n; kk = cos(2.0*pi*(k-1)/n)
        do j = 1, n; jj = cos(2.0*pi*(j-1)/n)
        do i = 1,nn; ii = cos(2.0*pi*(i-1)/n)
                Fk(i,j,k) = c*Fk(i,j,k)/(c0 + c1*(ii+jj+kk) + c2*(ii*jj+ii*kk+jj*kk) + c3*ii*jj*kk)
        end do; end do; end do
        
        Fk(1,1,1) = 0.0
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine laplace

! one-sided power spectrum density estimator
subroutine spectrum(f, S)
        integer, parameter :: ns = 1.732050807568877293527446341506*(n/2) + 1
        real f(n,n,n), S(ns), W(ns); complex Fk(nn,n,n); integer*8 plan
        integer i, j, k, ii, jj, kk, l; real p, c(2)
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,f,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        S = 0.0; W = 0.0
        
        do k = 1,n; if (k <= nn) then; kk = k-1; else; kk = n+1-k; end if
        do j = 1,n; if (j <= nn) then; jj = j-1; else; jj = n+1-j; end if
        do i = 1,n; if (i <= nn) then; ii = i-1; else; ii = n+1-i; end if
                p = sqrt(real(ii**2 + jj**2 + kk**2)); l = floor(p)
                
                !c = (/ 1.0, 0.0 /)
                c = (1.0 - (/l-p,l+1-p/)**2)**2
                
                S(l+1:l+2) = S(l+1:l+2) + c * Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
                W(l+1:l+2) = W(l+1:l+2) + c
        end do; end do; end do
        
        where (W .ne. 0.0) S = S/W/real(n)**6
end subroutine spectrum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
