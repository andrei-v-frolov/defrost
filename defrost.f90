! $Id: defrost.f90,v 1.15 2007/07/30 21:41:16 frolov Exp $
! [compile with: ifort -O3 -ipo -xT -r8 -pc80 -fpp defrost.f90 -lfftw3]

! Reheating code doing something...
! http://www.sfu.ca/physics/cosmology/defrost

! Copyright (C) 2007 Andrei Frolov <frolov@sfu.ca>
! Distributed under the terms of the GNU General Public License
! If you use this code for your research, please cite hep-../...

program defrost; implicit none

include "fftw3.f"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! some useful constants
real, parameter :: twopi = 6.2831853071795864769252867665590
real, parameter :: sqrt3 = 1.7320508075688772935274463415059

! solver control parameters
integer, parameter :: n = 32                   ! sampled grid size (simulation cube is n^3 pts)
integer, parameter :: p = n+2                   ! padded grid size (>n, adjust for cache efficiency)
integer, parameter :: tt = 2**18                ! total number of time steps to take (i.e. runtime)
integer, parameter :: nn = n/2+1                ! Nyquist frequency (calculated, leave it alone)
integer, parameter :: ns = sqrt3*(n/2) + 2      ! highest wavenumber on 3D grid (leave it alone)

real, parameter :: alpha = 40.0                 ! dx/dt (be careful not to violate Courant condition)
real, parameter :: dx = 1.0/n                  ! grid spacing   (physical grid size is n*dx)
real, parameter :: dt = dx/alpha                ! time step size (simulated timespan is tt*dt)
real, parameter :: dk = twopi/(n*dx)            ! frequency domain grid spacing (leave it alone)

! output control parameters
integer, parameter :: nx = n/2                  ! spatial grid is downsampled to nx^3 pts for output
integer, parameter :: nt = 2**9                 ! simulation will be logged every nt time steps

logical, parameter :: output = .true.           ! set this to false to disable all file output at once
logical, parameter :: oscale = .true.           ! scale output variables to counter-act expansion

logical, parameter :: output$bov = .false.      ! output 3D data cube (storage-expensive)
logical, parameter :: output$psd = .true.       ! output power spectra (time-expensive)
logical, parameter :: output$cdf = .true.       ! output distributions (time-expensive)

logical, parameter :: output$fld = .true.       ! output scalar fields
logical, parameter :: output$set = .true.       ! output stress-energy tensor components
logical, parameter :: output$pot = .true.       ! output gravitatinal potential (expensive)

logical, parameter :: output$gnu = .true.       ! output curves in gnuplot format (sinle file)
logical, parameter :: output$vis = .false.      ! output curves in VisIt X-Y format (per frame)
logical, parameter :: output$crv = (output$psd .or. output$cdf) .and. (output$gnu .or. output$vis)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fields are referred to by their symbolic aliases
integer, parameter :: fields = 2                ! total number of scalar fields being evolved
integer, parameter :: phi = 1, psi = 2          ! symbolic aliases for scalar fields
integer, parameter :: rho = 1, prs = 2          ! symbolic aliases for stress-energy

! potential and its derivatives are (separately) inlined in step()

! ... these will move ...
real, parameter :: m2phi = 1.0, m2psi = 0.0, g2 = 100.0**2, mpl = 2.0e5

! initial conditions for homogeneous field component
real, parameter ::  phi0 =  1.0093430384226378929425913902459
real, parameter :: dphi0 = -0.7137133070120812430962278466136
real, parameter ::    H0 =  0.5046715192113189464712956951230
real, parameter ::   dH0 = -H0**2
real, parameter :: ddphi0 = -(3.0*H0*dphi0 + m2phi*phi0)
real, parameter ::  ddH0 = -dphi0*ddphi0

! scale factor and horizon size (sampled on two subsequent time slices)
real :: LA(2) = (/ 1.0 - H0*dt, 1.0 /)
real :: LH(2) = 1.0/(/ H0 - dH0*dt + ddH0*dt**2/2.0, H0 /)

! ...
character(12) DVAR(fields+3)
real CDF(n+1,fields+3), PSD(ns,fields+3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer l

! smp array samples all fields on a 3D grid for three subsequent time slices
real smp(fields,0:p,0:p,0:p,3), tmp(n,n,n); complex Fk(nn,n,n)

! for larger grids, dynamically allocated array might be required
! real, allocatable :: smp(:,:,:,:,:), tmp(:,:,:); complex, allocatable :: Fk(:,:,:)
! allocate(smp(fields,0:p,0:p,0:p,3), tmp(n,n,n), Fk(nn,n,n))

! use threaded FFTW on SMP machines (link with -lfftw3_threads)
!call dfftw_init_threads
!call dfftw_plan_with_nthreads(4)

! initialize random number generator
call random_seed

! initialize and run simulation
call head(6, (/"t", "a", "H", "<rho>", "<P>"/))
call init(smp(:,:,:,:,1), smp(:,:,:,:,2))

do l = 1,tt,3
        call step(l,   smp(:,:,:,:,1), smp(:,:,:,:,2), smp(:,:,:,:,3), smp(:,:,:,:,1))
        call step(l+1, smp(:,:,:,:,2), smp(:,:,:,:,3), smp(:,:,:,:,1), smp(:,:,:,:,2))
        call step(l+2, smp(:,:,:,:,3), smp(:,:,:,:,1), smp(:,:,:,:,2), smp(:,:,:,:,3))
end do

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! periodic boundary conditions
subroutine wrap(up)
        real, dimension(fields,0:p,0:p,0:p) :: up
        
        up(:,0,:,:) = up(:,n,:,:); up(:,n+1,:,:) = up(:,1,:,:)
        up(:,:,0,:) = up(:,:,n,:); up(:,:,n+1,:) = up(:,:,1,:)
        up(:,:,:,0) = up(:,:,:,n); up(:,:,:,n+1) = up(:,:,:,1)
end subroutine wrap

! scalar field initial conditions
subroutine init(dn, hr)
        real, dimension(fields,0:p,0:p,0:p) :: dn, hr
        
        real, parameter :: m2phi$eff = m2phi - 2.25*H0**2
        real, parameter :: m2psi$eff = m2psi + g2*phi0**2 - 2.25*H0**2
        
        call sample(tmp, -0.25, m2phi$eff); hr(phi,1:n,1:n,1:n) = tmp + phi0
        call sample(tmp, -0.25, m2psi$eff); hr(psi,1:n,1:n,1:n) = tmp
        
        call sample(tmp, +0.25, m2phi$eff); dn(phi,1:n,1:n,1:n) = hr(phi,1:n,1:n,1:n) - tmp*dt - dphi0*dt + ddphi0*dt**2/2.0
        call sample(tmp, +0.25, m2psi$eff); dn(psi,1:n,1:n,1:n) = hr(psi,1:n,1:n,1:n) - tmp*dt
        
        call wrap(dn); call wrap(hr)
end subroutine init

! stencil operators (implemented as preprocessor macros)
#define HR(x,y,z) hr(:,i+(x),j+(y),k+(z))
#define GRAD2(x,y,z) sum((hr(:,i+(x),j+(y),k+(z))-hr(:,i,j,k))**2)

#define RANK0(O) (O(0,0,0))
#define RANK1(O) (O(-1,0,0) + O(1,0,0) + O(0,-1,0) + O(0,1,0) + O(0,0,-1) + O(0,0,1))
#define RANK2(O) (O(-1,-1,0) + O(1,-1,0) + O(-1,1,0) + O(1,1,0) + O(-1,0,-1) + O(1,0,-1) + O(-1,0,1) + O(1,0,1) + O(0,-1,-1) + O(0,1,-1) + O(0,-1,1) + O(0,1,1))
#define RANK3(O) (O(-1,-1,-1) + O(1,-1,-1) + O(-1,1,-1) + O(1,1,-1) + O(-1,-1,1) + O(1,-1,1) + O(-1,1,1) + O(1,1,1))
#define STENCIL(C,O) ((C ## 0) * RANK0(O) + (C ## 1) * RANK1(O) + (C ## 2) * RANK2(O) + (C ## 3) * RANK3(O))

! scalar field evolution step
subroutine step(l, dn, hr, up, pp)
        real, dimension(fields,0:p,0:p,0:p) :: dn, hr, up, pp; integer i, j, k, l
        
        ! Laplacian operator stencils: traditional one and three isotropic variants
        ! stable for dx/dt > sqrt(3), sqrt(2), sqrt(21)/3, and 8/sqrt(30) respectively
        ! computational cost difference is insignificant for large grids
        
        !real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -6.0, cc = 1.0
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 2.0, c0 = -24.0, cc = 6.0
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 8.0, c0 = -56.0, cc = 12.0
        real, parameter :: c3 = 1.0, c2 = 3.0, c1 = 14.0, c0 = -128.0, cc = 30.0
        
        ! basis vectors for vectorizing potential derivative
        real, dimension(fields), parameter :: V1 = (/1,0/), V2 = (/0,1/)
        
        ! field energy (distributed for parallelization)
        real, dimension(n) :: V, T, G, PE, KE, GE
        
        ! various evolution operator coefficients
        real c, d, b0, b1, b2, b3, d1, d2, e1, e2, e3
        
        ! optional computations flags
        logical, parameter :: dumping = output .and. (output$bov .or. output$crv)
        logical, parameter :: needTii = dumping .and. (output$set .or. output$pot)
        logical checkpt; integer idx
        
        ! flat or expanding background
        !real, parameter :: a = 1.0, H = 0.0; real Q
        real a, H, Q, R; a = LA(2); H = 1.0/LH(2)
        
        ! all coefficients inside the loop are pre-calculated here
        d = 1.0 + 1.5*H*dt; c = cc * alpha**2 * a**2 * d
        b0 = 2.0/d + c0/c; b1 = c1/c; b2 = c2/c; b3 = c3/c
        d1 = -(1.0 - 1.5*H*dt)/d; d2 = -dt**2/d
        e1 = 1.0/(8.0*dt**2); e2 = 1.0/(4.0*(a*dx)**2*cc); e3 = e2/3.0
        
        PE = 0.0; KE = 0.0; GE = 0.0
        checkpt = mod(l-1, nt) == 0; idx = 0
        
        do k = 1,n; do j = 1,n; do i = 1,n
                ! discretized scalar field evolution step
                ! scalar field potential derivatives are inlined here
                up(:,i,j,k) = STENCIL(b,HR) + d1 * dn(:,i,j,k) + &
                    d2 * ( (m2phi + g2*hr(psi,i,j,k)**2)*V1 + (m2psi + g2*hr(phi,i,j,k)**2)*V2 ) * hr(:,i,j,k)
                
                ! scalar field potential multiplied by 2 is inlined here
                V(k) = (m2phi + g2*hr(psi,i,j,k)**2)*hr(phi,i,j,k)**2 + m2psi*hr(psi,i,j,k)**2
                T(k) = sum((up(:,i,j,k)-dn(:,i,j,k))**2)
                
                PE(k) = PE(k) + V(k); KE(k) = KE(k) + T(k)
                
                ! calculate density and pressure when needed
                if (checkpt) then
                        G(k) = STENCIL(c,GRAD2); GE(k) = GE(k) + G(k)
                        
                        if (needTii) then
                                pp(rho,i,j,k) = e1*T(k) + e2*G(k) + 0.5*V(k)
                                pp(prs,i,j,k) = e1*T(k) - e3*G(k) - 0.5*V(k)
                        end if
                end if
        end do; end do; end do
        
        ! periodic boundary conditions
        call wrap(up)
        
        ! update expansion factors
        Q = sum(4.0*e1*KE - PE)/(6.0*n**3)
        R = LH(1) + (1.0 + Q * LH(2)**2) * dt
        LH = (/ LH(2), LH(1) + (1.0 + Q * R**2) * (2.0*dt) /)
        LA = (/ LA(2), LA(1) + (H*a) * (2.0*dt) /)
        
        ! dump simulation data
        if (checkpt) then
                write (*,'(5g)') (l-1)*dt, a, H, sum(e1*KE + e2*GE + 0.5*PE)/n**3, sum(e1*KE - e3*GE - 0.5*PE)/n**3
                
                if (dumping .and. output$fld) then
                        Q = 1.0; if (oscale) Q = a**1.5
                        tmp = Q*hr(phi,1:n,1:n,1:n); call dump("phi", (l-1)/nt, (l-1)*dt, tmp, idx)
                        tmp = Q*hr(psi,1:n,1:n,1:n); call dump("psi", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (dumping .and. output$set) then
                        Q = 1.0; if (oscale) Q = 1.0/(3.0*H**2)
                        tmp = Q*pp(rho,1:n,1:n,1:n); call dump("rho", (l-1)/nt, (l-1)*dt, tmp, idx)
                        tmp = Q*pp(prs,1:n,1:n,1:n); call dump("prs", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (dumping .and. output$pot) then
                        Q = a**2/2.0; tmp = Q*pp(rho,1:n,1:n,1:n)
                        call laplace(tmp, tmp); call dump("PSI", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (idx > 0 .and. output$gnu) call fflush((l-1)*dt, idx)
        end if
end subroutine step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sample Gaussian random field with power-law spectrum
subroutine sample(f, gamma, m2eff)
        real f(n,n,n), gamma, m2eff; integer*8 plan
        
        integer, parameter :: os = 16, nos = n * os**2
        real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = nn*dk/2.0
        real, parameter :: norm = 0.5/(n**3 * (twopi*dk**3)**0.5 * mpl) * (dkos/dxos)
        complex, parameter :: w = (0.0, twopi)
        
        real ker(nos), a(nn), p(nn)
        integer i, j, k, l; real kk
        
        ! calculate (oversampled) radial profile of convolution kernel
        do k = 1,nos; kk = (k-0.5)*dkos
                ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/kcut)**2)
        end do
        
        call dfftw_plan_r2r_1d(plan,nos,ker,ker,FFTW_RODFT10,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        do k = 1,nos; ker(k) = norm * ker(k)/k; end do
        
        ! initialize 3D convolution kernel (using linear interpolation of radial profile)
        do k = 1,n; do j = 1,n; do i = 1,n
                kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os; l = floor(kk)
                
                if (l > 0) then
                        f(i,j,k) = ker(l) + (kk-l)*(ker(l+1)-ker(l))
                else
                        f(i,j,k) = (4.0*ker(1)-ker(2))/3.0
                end if
        end do; end do; end do
        
        ! convolve kernel with delta-correlated Gaussian noise
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,f,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        do k = 1,n; do j = 1,n
                call random_number(a); call random_number(p)
                Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
        end do; end do
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine sample

! solve Laplace equation $\Delta f = \rho$
subroutine laplace(f, rho)
        real f(n,n,n), rho(n,n,n); integer*8 plan
        integer i, j, k; real :: ii, jj, kk
        
        real, parameter :: w = twopi/n
        
        ! Laplacian operator stencils: traditional one and three isotropic variants
        ! (corresponding to discretization used for field evolution equations above)
        
        !real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -3.0, cc = 0.5
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 1.0, c0 = -6.0, cc = 1.5
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 2.0, c0 = -7.0, cc = 1.5
        real, parameter :: c3 = 2.0, c2 = 3.0, c1 = 7.0, c0 = -32.0, cc = 7.5
        
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
        real f(n,n,n), S(ns), W(ns); integer*8 plan
        integer i, j, k, ii, jj, kk, l; real p, c(2)
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,f,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        S = 0.0; W = 0.0
        
        do k = 1,n; if (k <= nn) then; kk = k-1; else; kk = n+1-k; end if
        do j = 1,n; if (j <= nn) then; jj = j-1; else; jj = n+1-j; end if
        do i = 1,n; if (i <= nn) then; ii = i-1; else; ii = n+1-i; end if
                p = sqrt(real(ii**2 + jj**2 + kk**2)); l = floor(p)
                
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

! identify yourself
subroutine head(fd, vars)
        integer(4) fd; character(*) :: vars(:)
        character(512) :: buffer; integer a, b, l
        character(*), parameter :: rev = "$Revision: 1.15 $"
        
        a = index(rev, ": ") + 2
        b = index(rev, " $", .true.) - 1
        
        ! ID string and model summary
        write (fd,'(3g,2(i0,g))') "# This is DEFROST revision ", rev(a:b), " (", fields, " fields, ", n, "^3 grid)"
        write (fd,'(g,3(f0.5,g))') "# V(phi,psi) = ", &
                sqrt(m2phi), "^2*phi^2/2 + ",         &
                sqrt(m2psi), "^2*psi^2/2 + ",         &
                sqrt(g2), "^2*phi^2*psi^2/2"
        
        ! variable list
        write (buffer,'(g,g12.12",",8(g24.12","))') "OUTPUT:", adjustr(vars); l = index(buffer, ',', .true.);
        write (fd,'(2g)') "# ", repeat('=', l)
        write (fd,'(2g)') "# ", buffer(1:l-1)//';'
        write (fd,'(2g)') "# ", repeat('=', l)
end subroutine head

! output field configuration and its aggregates
subroutine dump(file, frame, t, f, idx)
        character(*) :: file; integer frame, k; real t, f(n,n,n); integer, optional :: idx
        character(256) :: buffer; real avg, var, S(ns), P(n+1), X(n+1);
        
        integer, parameter :: stride = n/nx
        
        ! output 3D box of values
        if (output$bov) then
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
                if (stride > 1) then; do k=1,n,stride; write (11) f(::stride,::stride,k); end do; else; write (11) f; end if
                
                close (10); close (11)
        end if
        
        ! output spectra and statistics
        if (output$crv) then
                ! in VisIt format, all curves go into a single file per field, per frame
                if (output$vis) then
                        write (buffer,'(a,a,i6.6,a)') file, '-', frame, '.ult'; open(12, file=buffer)
                end if
                
                ! in gnuplot format, all fields and frames go into a single file per curve
                ! (to be fflush()ed every frame after all the fields are analyzed)
                if (present(idx)) then; idx = idx + 1; DVAR(idx) = file; end if
                
                ! output power spectrum
                if (output$psd) then
                        call spectrum(f, S); if (present(idx)) PSD(:,idx) = S
                        
                        if (output$vis) then
                                write (12,'(g)') "# PSD"
                                do k = 1,ns; write (12,'(2g)') (k-1)*dk, S(k); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PSD [logarithmic]"
                                do k = 2,ns; write (12,'(2g)') log10((k-1)*dk), log10(S(k)); end do
                                write (12,'(g)') "", ""
                        end if
                end if
                
                ! output distribution of values
                if (output$cdf) then
                        call sieve(f, n**3, n**2, 1)
                        P(1:n) = f(1,1,:); P(n+1) = maxval(f(:,:,n)); if (present(idx)) CDF(:,idx) = P
                        avg = sum(P)/(n+1); var = sum((P-avg)**2)/n; X = (P-avg)/sqrt(2.0*var)
                        
                        if (output$vis) then
                                write (12,'(g)') "# CDF"
                                do k = 1,n+1; write (12,'(2g)') P(k), real(k-1)/n; end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# CDF [gaussian]"
                                do k = 1,n+1; write (12,'(2g)') P(k), (1.0 + erf(X(k)))/2.0; end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PDF"
                                do k = 2,n; write (12,'(2g)') P(k), (2.0/n)/(P(k+1)-P(k-1)); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PDF [gaussian]"
                                do k = 2,n; write (12,'(2g)') P(k), exp(-X(k)**2)/sqrt(twopi*var); end do
                                write (12,'(g)') "", ""
                        end if
                end if
                
                if (output$vis) close (12)
        end if
end subroutine dump

! flush frame data into gnuplot-style curves
subroutine fflush(t, idx)
        real t; integer idx, k; logical o
        
        ! output power spectrum
        if (output$psd) then
                inquire (31, opened=o)
                
                if (.not. o) then
                        open(31, file="PSD"); call head(31, (/"t", "k", DVAR(1:idx)/))
                end if
                
                do k = 1,ns; write (31,'(12g)') t, (k-1)*dk, log10(PSD(k,1:idx)); end do
                write (31,'(g)') "", ""; flush(31)
        end if
        
        ! output distribution of values
        if (output$cdf) then
                inquire (32, opened=o)
                
                if (.not. o) then
                        open(32, file="CDF"); call head(32, (/"t", "percentile", DVAR(1:idx)/))
                end if
                
                do k = 1,n+1; write (32,'(12g)') t, real(k-1)/n, CDF(k,1:idx); end do
                write (32,'(g)') "", ""; flush(32)
        end if
end subroutine fflush

end
