! $Id: defrost.f90,v 1.3 2007/05/05 06:51:48 frolov Exp $
! [compile with: ifc -O3 -ipo -xT -r8 -pc80 defrost.f90 -lfftw3]

! Reheating code doing something...

program main; implicit none

include "fftw3.f"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! some useful constants
real, parameter :: twopi = 6.2831853071795864769252867665590
real, parameter :: sqrt3 = 1.7320508075688772935274463415059

! control parameters
integer, parameter :: n = 256                   ! grid size: spatial grid contains n^3 points
integer, parameter :: p = n+2                   ! padded grid size (adjust for cache efficiency)
integer, parameter :: nn = n/2+1                ! Nyquist frequency
integer, parameter :: ns = sqrt3*(n/2) + 2      ! ...

integer, parameter :: nt = 32                   ! ...
integer, parameter :: nx = 32                   ! ...

real, parameter :: alpha = 2.0
real, parameter :: dx = 1.0/n, dt = dx/alpha


integer l

! smp array samples all fields on a 3D grid for three subsequent time slices
integer, parameter :: phi = 1, psi = 2, fields = 2; real smp(fields,0:p,0:p,0:p,3), tmp(n,n,n)

! for larger grids, dynamically allocated array might be required
! real, allocatable :: smp(:,:,:,:,:), tmp(:,:,:); allocate(smp(fields,0:p,0:p,0:p,3), tmp(n,n,n))

smp = 0.0
smp(1,n/2,n/2,n/2,:) = 1.0

!tmp = smp(1,1:n,1:n,1:n,1)
call sample(tmp, 1.0)
call dump("pHi", 1, (1-1)*dt, tmp)

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
        ! if (mod(l-1, n/nt) == 0) then; tmp = hr(phi,1:n,1:n,1:n); call dump("pHi", (l-1)/(n/nt), (l-1)*dt, tmp); end if
end subroutine step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sample Gaussian random field
subroutine sample(f, sigma)
        real f(n,n,n), sigma; complex Fk(nn,n,n); integer*8 plan
        integer i, j, k; real :: s = 0.0, a(nn), p(nn)
        
        complex, parameter :: w = (0.0, twopi)
        
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
        
        real, parameter :: w = twopi/n
        
        real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -3.0, cc = 0.5
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 1.0, c0 = -6.0, cc = 1.5
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 2.0, c0 = -7.0, cc = 1.5
        !real, parameter :: c3 = 2.0, c2 = 3.0, c1 = 7.0, c0 = -32.0, cc = 7.5
        
        real, parameter :: c = cc * dx**2/real(n)**3
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,rho,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        do k = 1, n; kk = cos(w*(k-1))
        do j = 1, n; jj = cos(w*(j-1))
        do i = 1,nn; ii = cos(w*(i-1))
                Fk(i,j,k) = c*Fk(i,j,k)/(c0 + c1*(ii+jj+kk) + c2*(ii*jj+ii*kk+jj*kk) + c3*ii*jj*kk)
        end do; end do; end do
        
        Fk(1,1,1) = 0.0
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine laplace

! one-sided power spectrum density estimator
subroutine spectrum(f, S)
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
        
        where (W /= 0.0) S = S/W/real(n)**6
end subroutine spectrum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! partially sort an array, finding (r:n:m)th smallest elements
recursive subroutine sieve(f, n, m, r)
        integer i, j, k, n, m, r; real p, f(n)
        
        if (r > n) return
        
        i = 1; k = n/2+1; j = n
        
        if (f(1) > f(n)) f((/1,n/)) = f((/n,1/))
        if (f(1) > f(k)) f((/1,k/)) = f((/k,1/))
        if (f(n) < f(k)) f((/n,k/)) = f((/k,n/))
        
        if (n < 4) return
        
        p = f(k)
        
        do while (i < j)
                i = i+1; do while (f(i) < p); i = i+1; end do
                j = j-1; do while (f(j) > p); j = j-1; end do
                if (i < j) f((/i,j/)) = f((/j,i/))
        end do
        
        call sieve(f(1:i-1), i-1, m, r)
        call sieve(f(j+1:n), n-j, m, modulo(r-j-1,m)+1)
end subroutine sieve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 
subroutine dump(file, frame, t, f)
        character(*) :: file; integer frame, k; real t, f(n,n,n);
        character(256) :: buffer; real avg, var, S(ns), P(n+1), X(n+1);
        
        integer, parameter :: stride = n/nx
        real, parameter :: dk = twopi/(n*dx)
        
        write (buffer,'(a,a,i6.6,a)') file, '-', frame, '.ult'; open(12, file=buffer)
        write (buffer,'(a,a,i6.6,a)') file, '-', frame, '.bov'; open(10, file=buffer)
        write (buffer,'(a,a,i6.6,a)') file, '-', frame, '.raw'; open(11, file=buffer, form="binary")
        
        ! BOV header
        write (10,'(4g)') "TIME:        ", t
        write (10,'(4g)') "BRICK_ORIGIN:", 0.0, 0.0, 0.0
        write (10,'(4g)') "BRICK_SIZE:  ", n*dx, n*dx, n*dx
        write (10,'(4g)') "VARIABLE:        ", file
        write (10,'(4g)') "DATA_FILE:       ", buffer
        write (10,'(g,3i5)') "DATA_SIZE:     ", nx, nx, nx
        write (10,'(4g)') "DATA_FORMAT:     ", "DOUBLE"
        write (10,'(4g)') "DATA_ENDIAN:     ", "LITTLE"
        write (10,'(4g)') "CENTERING:       ", "zonal"
        
        ! raw data
        if (stride > 1) then; write (11) f(::stride,::stride,::stride); else; write (11) f; end if
        
        ! power spectrum
        call spectrum(f, S)
        
        write (12,'(g)') "# PSD"
        do k = 1,ns; write (12,'(2g)') (k-1)*dk, S(k); end do
        
        write (12,'(g)') "", "", "# PSD [logarithmic]"
        do k = 2,ns; write (12,'(2g)') log10((k-1)*dk), log10(S(k)); end do
        
        ! value distribution
        call sieve(f, n**3, n**2, 1)
        P(1:n) = f(1,1,:); P(n+1) = maxval(f(:,:,n))
        avg = sum(P)/(n+1); var = sum((P-avg)**2)/n; X = (P-avg)/sqrt(2.0*var)
        
        write (12,'(g)') "", "", "# CDF"
        do k = 1,n+1; write (12,'(2g)') P(k), real(k-1)/n; end do
        
        write (12,'(g)') "", "", "# CDF [gaussian]"
        do k = 1,n+1; write (12,'(2g)') P(k), (1.0 + erf(X(k)))/2.0; end do
        
        write (12,'(g)') "", "", "# PDF"
        do k = 2,n; write (12,'(2g)') P(k), (2.0/n)/(P(k+1)-P(k-1)); end do
        
        write (12,'(g)') "", "", "# PDF [gaussian]"
        do k = 2,n; write (12,'(2g)') P(k), exp(-X(k)**2)/sqrt(twopi*var); end do
        
        close (10); close (11); close (12)
end subroutine dump

end
