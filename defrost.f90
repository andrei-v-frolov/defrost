! $Id: defrost.f90,v 3.0 2012/04/21 09:36:43 frolov Exp $
! [compile with: ifort -O3 -ipo -xT -r8 -pc80 -fpp defrost.f90 -lfftw3]

! Reheating code doing something...
! http://www.sfu.ca/physics/cosmology/defrost

! Copyright (C) 2007-2012 Andrei Frolov <frolov@sfu.ca>
! Distributed under the terms of the GNU General Public License
! If you use this code for your research, please cite arXiv:0809.4904

program defrost; implicit none

include "fftw3.f"

#ifdef SILO
include "silo.inc"
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preheating model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fields are referred to by their symbolic aliases
integer, parameter :: fields = 2                ! total number of scalar fields being evolved
integer, parameter :: phi = 1, psi = 2          ! symbolic aliases for scalar field components

! parameters of the scalar field potential are defined here
real, parameter :: lambda = 1.0, g2 = 1.875, mpl = 1.0e7/3.0

! potential and its derivatives are inlined in ...step() routines
#define Vx4(PHI,PSI) (lambda * (PHI)**2 + (2.0*g2) * (PSI)**2) * (PHI)**2
#define M2I(PHI,PSI) VECTOR(lambda*(PHI)**2 + g2*(PSI)**2, g2*(PHI)**2)

! model summary in human-readable form is printed out in head()
#define VTAG "# V(phi,psi) = ", lambda, "/4 * phi^4 + ", g2, "/2 phi^2 psi^2"

! initial conditions for homogeneous field components
real, parameter ::  phi0 =  2.339383796213256
real, parameter :: dphi0 = -2.736358272992573
real, parameter ::    H0 =  1.934897490588959


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Configurable parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! some useful constants (in quad precision)
real, parameter :: twopi = 6.28318530717958647692528676655900577
real, parameter :: sqrt3 = 1.73205080756887729352744634150587237

! solver control parameters
integer, parameter :: n = 256                   ! sampled grid size (simulation cube is n^3 pts)
integer, parameter :: m = n+2                   ! padded grid size (>n, adjust for cache efficiency)
integer, parameter :: tt = 2**13                ! total number of time steps to take (i.e. runtime)
integer, parameter :: nn = n/2+1                ! Nyquist frequency (calculated, leave it alone)
integer, parameter :: ns = sqrt3*(n/2) + 2      ! highest wavenumber on 3D grid (leave it alone)

real, parameter :: alpha = 2.5                  ! dx/dt (be careful not to violate Courant condition)
real, parameter :: dx = 10.0/n                  ! grid spacing   (physical grid size is n*dx)
real, parameter :: dt = dx/alpha                ! time step size (simulated timespan is tt*dt)
real, parameter :: dk = twopi/(n*dx)            ! frequency domain grid spacing (leave it alone)

! output control parameters
integer, parameter :: nx = n/2                  ! spatial grid is downsampled to nx^3 pts for output
integer, parameter :: nt = 2**4                 ! simulation will be logged every nt time steps

logical, parameter :: output = .true.           ! set this to false to disable all file output at once
logical, parameter :: oscale = .true.           ! scale output variables to counter-act expansion

logical, parameter :: output$bov = .false.      ! output 3D data cube (storage-expensive)
logical, parameter :: output$psd = .true.       ! output power spectra (time-expensive)
logical, parameter :: output$cdf = .true.       ! output distributions (time-expensive)

logical, parameter :: output$fld = .true.       ! output scalar field values
logical, parameter :: output$vel = .false.      ! output scalar field velocities
logical, parameter :: output$rho = .true.       ! output total energy density
logical, parameter :: output$prs = .true.       ! output isotropic pressure
logical, parameter :: output$pot = .true.       ! output gravitatinal potential (expensive)

logical, parameter :: output$gnu = .true.       ! output curves in gnuplot format (sinle file)
logical, parameter :: output$vis = .false.      ! output curves in VisIt X-Y format (per frame)
logical, parameter :: output$crv = (output$psd .or. output$cdf) .and. (output$gnu .or. output$vis)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Numerical scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! time integration is now symplectic, scheme order can be adjusted
! unless you have specific reasons to, stick to the 6-th order one
integer, parameter :: order = 6 ! integrator order (should be even)

! higher order symplectic integrator schedules (these are magical, don't touch)
real, parameter :: W6A(0:3) = (/ 1.31518632068391121888424972823886251, -1.17767998417887100694641568096431573, &
                                 0.235573213359358133684793182978534602, 0.784513610477557263819497633866349876 /)

! Laplacian operator stencils: traditional one and three isotropic variants
! stable for dx/dt > sqrt(3), sqrt(2), sqrt(21)/3, and 8/sqrt(30) respectively
! computational cost difference is insignificant for large grids

!character, parameter :: scheme = 'N'; real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -6.0, cc = 1.0
!character, parameter :: scheme = 'A'; real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 2.0, c0 = -24.0, cc = 6.0
!character, parameter :: scheme = 'B'; real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 8.0, c0 = -56.0, cc = 12.0
character, parameter :: scheme = 'C'; real, parameter :: c3 = 1.0, c2 = 3.0, c1 = 14.0, c0 = -128.0, cc = 30.0

! stencil operators (implemented as preprocessor macros)
#define RANK0(O) (O(0,0,0))
#define RANK1(O) (O(-1,0,0) + O(1,0,0) + O(0,-1,0) + O(0,1,0) + O(0,0,-1) + O(0,0,1))
#define RANK2(O) (O(-1,-1,0) + O(1,-1,0) + O(-1,1,0) + O(1,1,0) + O(-1,0,-1) + O(1,0,-1) + O(-1,0,1) + O(1,0,1) + O(0,-1,-1) + O(0,1,-1) + O(0,-1,1) + O(0,1,1))
#define RANK3(O) (O(-1,-1,-1) + O(1,-1,-1) + O(-1,1,-1) + O(1,1,-1) + O(-1,-1,1) + O(1,-1,1) + O(-1,1,1) + O(1,1,1))
#define STENCIL(C,O) ((C ## 0) * RANK0(O) + (C ## 1) * RANK1(O) + (C ## 2) * RANK2(O) + (C ## 3) * RANK3(O))

! parallel-safe vector constructor (... explain this ...)
#ifndef THREADS
#define VECTOR(A,B) (/ A, B /)
#else
#define VECTOR(A,B) ((A)*(/1,0/) + (B)*(/0,1/))
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer l, seed(2)

! canonical expansion variables
real :: q = 1.0, p = -6.0*H0

! buffers holding variable names, statistics and spectra (on per-frame basis)
character(12) DVAR(fields*2+3); real CDF(n+1,fields*2+3), PSD(ns,fields*2+3)

! allocate global storage for 3D grid variables and FFTs
! for larger grids, dynamically allocated array might be required
#ifndef DYNAMIC_ARRAYS
real fi(fields,0:m,0:m,0:m), pi(fields,0:m,0:m,0:m), rho(n,n,n), prs(n,n,n)
real(8) tmp(n,n,n), bov(nx,nx,nx); complex(8) Fk(nn,n,n)
#else
real, allocatable :: fi(:,:,:,:), pi(:,:,:,:), rho(:,:,:), prs(:,:,:)
real(8), allocatable :: tmp(:,:,:), bov(:,:,:)
complex(8), allocatable :: Fk(:,:,:)

allocate(fi(fields,0:m,0:m,0:m), pi(fields,0:m,0:m,0:m))
allocate(tmp(n,n,n), Fk(nn,n,n))
if (output$rho) allocate(rho(n,n,n))
if (output$prs) allocate(prs(n,n,n))
if (output$bov .and. nx /= n) allocate(bov(nx,nx,nx))
#endif

! use threaded FFTW on SMP machines (link with -lfftw3_threads)
#ifdef THREADS
call dfftw_init_threads(l)
call dfftw_plan_with_nthreads(THREADS)
#endif

! initialize random number generator (use urandom on clusters!)
open (333, file="/dev/random", action='read', form='binary')
read (333) seed; call random_seed(PUT=seed); close (333)

! initialize simulation
call head(6, (/"t", "a", "H", "<rho>", "<P>", "Omega[T]", "Omega[G]", "Omega[V]", "Omega[K]", "<phi>", "<psi>"/))
call init(fi, pi); call checkpt(0, fi, pi)

! time evolution loop
do l = 1,tt
        call si(fi, pi, dt, order); if (mod(l, nt) == 0) call checkpt(l, fi, pi)
end do

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! periodic boundary conditions
subroutine wrap(up)
        real, dimension(fields,0:m,0:m,0:m) :: up
        
        up(:,0,:,:) = up(:,n,:,:); up(:,n+1,:,:) = up(:,1,:,:)
        up(:,:,0,:) = up(:,:,n,:); up(:,:,n+1,:) = up(:,:,1,:)
        up(:,:,:,0) = up(:,:,:,n); up(:,:,:,n+1) = up(:,:,:,1)
end subroutine wrap

! scalar field initial conditions
subroutine init(fi, pi)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi
        
        ! effective mass of the field fluctuations
        real, parameter :: m2phi$eff = 3.0*lambda*phi0**2
        real, parameter :: m2psi$eff = g2*phi0**2
        
        ! initialize simulation volume
        call sample(1, tmp, m2phi$eff, H0, phi0); fi(phi,1:n,1:n,1:n) = tmp
        call sample(1, tmp, m2psi$eff, H0      ); fi(psi,1:n,1:n,1:n) = tmp
        
        call sample(2, tmp, m2phi$eff, H0, dphi0); pi(phi,1:n,1:n,1:n) = tmp
        call sample(2, tmp, m2psi$eff, H0       ); pi(psi,1:n,1:n,1:n) = tmp
end subroutine init

! split Hamiltonian evolution - field advance step
subroutine fstep(fi, pi, dt, fuse)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real dt
        real, save :: d = 0.0; logical, optional :: fuse
        
        ! accumulate step
        d = d + dt/q**2
        
        ! advance fields (if needed)
        if (present(fuse) .and. fuse) return
        if (d /= 0.0) fi = fi + d*pi
        call wrap(fi); d = 0.0
end subroutine fstep

! split Hamiltonian evolution - momenta advance step
subroutine pstep(fi, pi, dt)
#define FI(x,y,z) fi(:,i+(x),j+(y),k+(z))
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real dt; integer i, j, k
        real, dimension(n) :: T, G, V; real, dimension(fields, n) :: DD, DV
        real a, b, b0, b1, b2, b3, d
        
        ! all coefficients inside the loop are pre-calculated here
        a = q; b = dt/cc * (a/dx)**2; d = a**4 * dt
        b0 = b*c0; b1 = b*c1; b2 = b*c2; b3 = b*c3;
        
        ! initialize accumulators
        T = 0.0; G = 0.0; V = 0.0
        
        !$omp parallel do
        do k = 1,n; do j = 1,n; do i = 1,n
                ! laplacian operator and scalar field potential derivatives
                DD(:,k) = STENCIL(b,FI)
                DV(:,k) = M2I(fi(phi,i,j,k),fi(psi,i,j,k)) * fi(:,i,j,k)
                
                ! accumulate field energy (pass 1)
                T(k) = T(k) + sum(pi(:,i,j,k)**2)
                G(k) = G(k) + sum(fi(:,i,j,k)*DD(:,k))
                
                ! advance field momenta
                pi(:,i,j,k) = pi(:,i,j,k) + DD(:,k) - d * DV(:,k)
                
                ! accumulate field energy (pass 2)
                V(k) = V(k) + Vx4(fi(phi,i,j,k),fi(psi,i,j,k))
                T(k) = T(k) + sum(pi(:,i,j,k)**2)
        end do; end do; end do
        
        ! reduce accumulated values and advance expansion rate
        p = p + sum(T)/(a*n)**3 * (dt/2.0) + sum(G)/(a*n**3) - sum(V)*(a/n)**3 * dt
end subroutine pstep

! split Hamiltonian evolution - logging step
subroutine lstep(fi, pi, time, flatten, need$rho, need$prs)
#define G2(x,y,z) (fi(:,i+(x),j+(y),k+(z))-fi(:,i,j,k))**2
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real time
        logical flatten, need$rho, need$prs; integer i, j, k
        real, dimension(n) :: T, G, V, KA, GA, PA
        real FA(fields,n), F(fields), KE, GE, PE, RE
        real a, H, e1, e2, e3, e4
        
        ! update fields to current time
        call fstep(fi, pi, 0.0); a = q
        
        ! all coefficients inside the loop are pre-calculated here
        e1 = 1.0/(2.0*a**4); e2 = 1.0/(4.0*cc*dx**2); e3 = e2/3.0; e4 = a**2/4.0
        
        ! initialize accumulators
        KA = 0.0; GA = 0.0; PA = 0.0; FA = 0.0
        
        !$omp parallel do
        do k = 1,n; do j = 1,n; do i = 1,n
                G(k) = sum(STENCIL(c,G2))
                T(k) = sum(pi(:,i,j,k)**2)
                V(k) = Vx4(fi(phi,i,j,k),fi(psi,i,j,k))
                
                FA(:,k) = FA(:,k) + fi(:,i,j,k)
                
                KA(k) = KA(k) + T(k)
                GA(k) = GA(k) + G(k)
                PA(k) = PA(k) + V(k)
                
                if (need$rho) rho(i,j,k) = e1*T(k) + e2*G(k) + e4*V(k)
                if (need$prs) prs(i,j,k) = e1*T(k) - e3*G(k) - e4*V(k)
        end do; end do; end do
        
        ! reduce accumulated values
        KE = e1*sum(KA)/n**3
        GE = e2*sum(GA)/n**3
        PE = e4*sum(PA)/n**3
        
        ! reset expansion rate to enforce flatness
        if (flatten) p = -sqrt(12.0*(KE+GE+PE))/a
        H = -p/(6.0*a); RE = KE+GE+PE-3.0*H**2
        
        forall (i=1:fields) F(i) = sum(FA(i,:))/n**3
        
        ! log expansion history and diagnostics to stdout
        write (*,'(32g)') time, a, H, KE+GE+PE, KE-GE/3.0-PE, (/KE, GE, PE, RE/)/(3.0*H**2), a*F
end subroutine lstep

! scalar field evolution checkpoint
subroutine checkpt(l, fi, pi)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; integer i, j, k, l, db, idx; real a, H, s
        
        ! optional computations flags
        logical, parameter :: dumping = output .and. (output$bov .or. output$crv) .and. &
                (output$fld .or. output$vel .or. output$rho .or. output$prs .or. output$pot)
        logical, parameter :: need$rho = dumping .and. (output$rho .or. output$pot)
        logical, parameter :: need$prs = dumping .and. output$prs
        
        ! prepare to log (flatten initial slice only!)
        call lstep(fi, pi, l*dt, l==0, need$rho, need$prs); a = q; H = -p/(6.0*a)
        
        ! dump requested variables
        db = 0; idx = 0
        if (dumping) db = fopen("frame", l/nt, l*dt)
        if (dumping .and. output$fld) then
                s = 1.0; if (oscale) s = a
                tmp = s*fi(phi,1:n,1:n,1:n); call dump(db, "phi", l/nt, l*dt, tmp, idx)
                tmp = s*fi(psi,1:n,1:n,1:n); call dump(db, "psi", l/nt, l*dt, tmp, idx)
        end if
        if (dumping .and. output$vel) then
                s = 1.0/a**2; if (oscale) s = 1.0/a
                tmp = s*pi(phi,1:n,1:n,1:n); call dump(db, "dphi", l/nt, l*dt, tmp, idx)
                tmp = s*pi(psi,1:n,1:n,1:n); call dump(db, "dpsi", l/nt, l*dt, tmp, idx)
        end if
        if (dumping .and. output$rho) then
                s = 1.0; if (oscale) s = 1.0/(3.0*H**2); tmp = s*rho
                call dump(db, "rho", l/nt, l*dt, tmp, idx)
        end if
        if (dumping .and. output$prs) then
                s = 1.0; if (oscale) s = 1.0/(3.0*H**2); tmp = s*prs
                call dump(db, "prs", l/nt, l*dt, tmp, idx)
        end if
        if (dumping .and. output$pot) then
                tmp = 0.5*rho; call laplace(tmp, tmp)
                call dump(db, "PSI", l/nt, l*dt, tmp, idx)
        end if
        if (idx > 0 .and. output$gnu) call fflush(l*dt, idx)
        if (db /= 0) call fclose(db)
end subroutine checkpt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! second order symplectic integrator
subroutine si2(fi, pi, dt)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real dt
        
        q = q - p * (dt/12.0)
        call fstep(fi, pi, dt/2.0)
        call pstep(fi, pi, dt    )
        call fstep(fi, pi, dt/2.0, fuse=.true.)
        q = q - p * (dt/12.0)
end subroutine si2

! 6-th order symplectic integrator of Yoshida (scheme A)
subroutine si6(fi, pi, dt)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real dt; integer i
        
        do i = -3,3; call si2(fi, pi, W6A(abs(i))*dt); end do
end subroutine si6

! k-th order recursive symplectic integrator (k should be even)
recursive subroutine si(fi, pi, dt, k)
        real, dimension(fields,0:m,0:m,0:m) :: fi, pi; real dt, gamma, w1, w0; integer k
        
        select case (k)
                case (2); call si2(fi, pi, dt)
                case (6); call si6(fi, pi, dt)
                case default
                        gamma = 1.0/(k-1); w1 = 1.0/(2.0 - 2.0**gamma); w0 = 1.0 - 2.0*w1
                        
                        call si(fi, pi, w1*dt, k-2)
                        call si(fi, pi, w0*dt, k-2)
                        call si(fi, pi, w1*dt, k-2)
        end select
end subroutine si


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! integrate modulus of a Hankel function w(x) = sqrt(pi/2x) |H_nu(x)|
! w(x) satisfies ODE w'' + (1+mu/x^2) w = 1/w^3, with mu = 1/4-nu^2
! see http://dlmf.nist.gov/10.18 for details and asymptotic expansion
function hankelw(y, mu, x, kind)
        real hankelw, y(3), mu, x, z, u, v, w, delta, dz; integer i, j, s, kind
        
        ! initialize based on asymptotic expansion at large x
        if (all(y == 0.0)) then
                z = 1.0/x
                w = sqrt(1.0 - (mu/2.0)*z**2)
                u = z*w; v = (1.0 - mu*z**2)/w
        else
                z = y(1); u = y(2); v = y(3)
        end if
        
        ! estimate number of steps needed to achive target accuracy
        delta = 1.0/x-z; s = 1 + int(abs(delta)/0.01); delta = delta/s
        
        ! take a sequence of leapfrog steps, 6th order schedule
        do j=1,s; do i = -3,3; dz = W6A(abs(i))/2.0 * delta
                z = z + dz; u = u + v * dz
                v = v + (1.0/u**3 - (1.0/z**2 + mu) * u/z**2) * (2.0*dz)
                z = z + dz; u = u + v * dz
        end do; end do
        
        ! field and velocity correction factors
        select case (kind)
                case(1); w = u/z
                case(2); w = sqrt((z*v - 2.0*u)**2 + (z/u)**2)
        end select
        
        ! pack dynamical system state and return
        y = (/z, u, v/); hankelw = w
end function hankelw

! sample Gaussian random field with specified spectrum
subroutine sample(kind, f, m2eff, H0, f0)
        real(8) f(n,n,n); integer(8) plan
        integer, value :: kind; real m2eff, H0; real, optional :: f0
        
        integer, parameter :: os = 16, nos = n * os**2
        real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = nn*dk/2.0
        real, parameter :: norm = 0.5/(n**3 * (twopi*dk**3)**0.5 * mpl) * (dkos/dxos)
        complex, parameter :: w = (0.0, twopi)
        
        real(8) ker(nos); real a(nn), p(nn)
        integer i, j, k, l; real kk, mu, y(3); y = 0.0
        
        ! vacuum state selection
        if (kind < 3 .and. H0 == 0.0) kind = kind + 2
        select case (kind)
                case (1:2); mu = m2eff/H0**2 - 2.0
                case (3:4); mu = m2eff - 2.25*H0**2
        end select
        
        ! initial fluctuations spectrum
        do k = nos,1,-1; kk = (k-0.5)*dkos
                select case (kind)
                        ! de Sitter vacuum
                        case(1); ker(k) = hankelw(y, mu, kk/H0, 1)/sqrt(kk)
                        case(2); ker(k) = hankelw(y, mu, kk/H0, 2)*sqrt(kk)
                        ! Minkowski vacuum
                        case(3); ker(k) = (mu + kk**2)**(-0.25)
                        case(4); ker(k) = (mu + kk**2)**(+0.25)
                end select
                
                ! cutoff is necessary for stability
                ker(k) = kk*ker(k) * exp(-(kk/kcut)**2)
        end do
        
        ! calculate (oversampled) radial profile of convolution kernel
        call dfftw_plan_r2r_1d(plan,nos,ker,ker,FFTW_RODFT10,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        do k = 1,nos; ker(k) = norm * ker(k)/k; end do
        
        ! initialize 3D convolution kernel (using linear interpolation of radial profile)
        !$omp parallel do
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
        
        !$omp parallel do
        do k = 1,n; do j = 1,n
                call random_number(a); call random_number(p)
                Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
        end do; end do
        
        ! override zero mode if requested
        if (present(f0)) Fk(1,1,1) = f0
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine sample

! solve Laplace equation $\Delta f = \rho$
subroutine laplace(f, rho)
        real(8) f(n,n,n), rho(n,n,n); integer(8) plan
        integer i, j, k; real :: ii, jj, kk
        
        real, parameter :: w = twopi/n, b = cc * dx**2/real(n)**3
        real, parameter :: b3 = 8.0*c3, b2 = 4.0*c2, b1 = 2.0*c1, b0 = c0
        
        call dfftw_plan_dft_r2c_3d(plan,n,n,n,rho,Fk,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        !$omp parallel do
        do k = 1, n; kk = cos(w*(k-1))
        do j = 1, n; jj = cos(w*(j-1))
        do i = 1,nn; ii = cos(w*(i-1))
                Fk(i,j,k) = b*Fk(i,j,k)/(b0 + b1*(ii+jj+kk) + b2*(ii*jj+ii*kk+jj*kk) + b3*ii*jj*kk)
        end do; end do; end do
        
        Fk(1,1,1) = 0.0
        
        call dfftw_plan_dft_c2r_3d(plan,n,n,n,Fk,f,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
end subroutine laplace

! one-sided power spectrum density estimator
subroutine spectrum(f, S)
        real(8) f(n,n,n), S(ns), W(ns); integer(8) plan
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

! partially sort an array, finding (r:n:m)th smallest elements
recursive subroutine sieve(f, n, m, r)
        integer i, j, k, n, m, r; real(8) p, f(n)
        
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
! I/O Backend
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! identify yourself
subroutine head(fd, vars)
        integer(4) fd; character(*) :: vars(:)
        character(512) :: buffer; integer a, b, c, l
        character(*), parameter :: rev = "$Revision: 3.0 $"
        
        ! behold the horror that is Fortran string parsing
        a = index(rev, ": ") + 2
        b = index(rev, ".")
        c = index(rev, " $", .true.) - 1
        buffer = rev(a:b-1); read (buffer, '(i)') a
        buffer = rev(b+1:c); read (buffer, '(i)') b
        
        ! ID string
        write (fd,'(g,5(i0,g))') "# This is DEFROST revision ", a-1, ".", b, " (", fields, " fields, ", n, "^3 grid, ", tt, " steps)"
        write (fd,'(g,i0,g,g)') "# ", order, "-th order symplectic integrator, Laplacian discretization ", scheme
        
        ! model summary
        write (fd,'(g,4(f0.5,g))') VTAG
        
        ! variable list
        write (buffer,'(g,g12.12",",32(g24.12","))') "OUTPUT:", adjustr(vars); l = index(buffer, ',', .true.);
        write (fd,'(2g)') "# ", repeat('=', l)
        write (fd,'(2g)') "# ", buffer(1:l-1)//';'
        write (fd,'(2g)') "# ", repeat('=', l)
end subroutine head

! open frame
function fopen(file, frame, t)
        character(*) :: file; integer fopen, frame; real t

#ifdef SILO
        integer db, opts, e
        integer, parameter :: D = 3, nodes(D) = nx+1
        integer, parameter :: lut(nx) = (/0:nx-1/)*n/nx
        real(8), parameter :: mesh(nx+1) = (/lut*dx,n*dx/)
        
        character(256) :: buffer; write (buffer,'(a,a,i4.4,a)') file, '-', frame, '.silo'
        
        if (.not. (output$bov .or. output$vis)) return
        
#define STR(string) string, len(string)
        e = dbcreate(buffer, len_trim(buffer), DB_CLOBBER, DB_LOCAL, STR("DEFROST frame"), DB_HDF5, db)
        
        e = dbmkoptlist(3, opts)
        e = dbaddiopt(opts, DBOPT_CYCLE, frame)
        e = dbadddopt(opts, DBOPT_DTIME, real(t,kind=8))
        e = dbaddiopt(opts, DBOPT_COORDSYS, DB_CARTESIAN)
        
        e = dbputqm(db, STR("mesh"), STR("x"), STR("y"), STR("z"), mesh, mesh, mesh, nodes, D, DB_DOUBLE, DB_COLLINEAR, opts, e)
        e = dbfreeoptlist(opts)
        
        fopen = db
#else
        fopen = 0
#endif
end function fopen

! close frame
subroutine fclose(db)
        integer db, e
        
#ifdef SILO
        e = dbclose(db)
#endif
end subroutine fclose

! output field configuration and its aggregates
subroutine dump(db, v, frame, t, f, idx)
        character(*) :: v; integer db, frame, k; real t; real(8) f(n,n,n); integer, optional :: idx
        
        character(256) :: buffer
        real(8) S(ns), C(n+1), PDF(n+1), X(n+1), Y(n+1), Z(n+1), avg, var
        real(8), parameter :: Q(ns) = (/0:ns-1/)*dk, P(n+1) = (/0:n/)/real(n)
        
        integer, parameter :: D = 3, zones(D) = nx, lut(nx) = (/0:nx-1/)*n/nx + 1; integer e
        
        ! output 3D box of values
        if (output$bov) then
#ifdef SILO
                if (nx /= n) then
                        bov = f(lut,lut,lut)
                        e = dbputqv1(db, STR(v), STR("mesh"), bov, zones, D, DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, e)
                else
                        e = dbputqv1(db, STR(v), STR("mesh"), f, zones, D, DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, e)
                end if
#else
                write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.bov'; open(10, file=buffer)
                write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.raw'; open(11, file=buffer, form="binary")
                
                ! BOV header
                write (10,'(4g)') "TIME:        ", t
                write (10,'(4g)') "BRICK_ORIGIN:", 0.0, 0.0, 0.0
                write (10,'(4g)') "BRICK_SIZE:  ", n*dx, n*dx, n*dx
                write (10,'(4g)') "VARIABLE:        ", v
                write (10,'(4g)') "DATA_FILE:       ", buffer
                write (10,'(g,3i5)') "DATA_SIZE:     ", nx, nx, nx
                write (10,'(4g)') "DATA_FORMAT:     ", "DOUBLE"
                write (10,'(4g)') "DATA_ENDIAN:     ", "LITTLE"
                write (10,'(4g)') "CENTERING:       ", "zonal"
                
                ! raw data
                if (nx /= n) then; bov = f(lut,lut,lut); write (11) bov; else; write (11) f; end if
                
                close (10); close (11)
#endif
        end if
        
        ! output spectra and statistics
        if (output$crv) then
#ifndef SILO
                ! in VisIt format, all curves go into a single file per field, per frame
                if (output$vis) then
                        write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.ult'; open(12, file=buffer)
                end if
#endif
                
                ! in gnuplot format, all fields and frames go into a single file per curve
                ! (to be fflush()ed every frame after all the fields are analyzed)
                if (present(idx)) then; idx = idx + 1; DVAR(idx) = v; end if
                
                ! output power spectrum
                if (output$psd) then
                        call spectrum(f, S); if (present(idx)) PSD(:,idx) = S
                        
                        if (output$vis) then
#ifdef SILO
                                e = dbputcurve(db, STR("PSD_"//v), Q, S, DB_DOUBLE, ns, DB_F77NULL, e)
                                e = dbputcurve(db, STR("log10_PSD_"//v), log10(Q(2:ns)), log10(S(2:ns)), DB_DOUBLE, ns-1, DB_F77NULL, e)
#else
                                write (12,'(g)') "# PSD"
                                do k = 1,ns; write (12,'(2g)') Q(k), S(k); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PSD [logarithmic]"
                                do k = 2,ns; write (12,'(2g)') log10(Q(k)), log10(S(k)); end do
                                write (12,'(g)') "", ""
#endif
                        end if
                end if
                
                ! output distribution of values
                if (output$cdf) then
                        call sieve(f, n**3, n**2, 1)
                        C(1:n) = f(1,1,:); C(n+1) = maxval(f(:,:,n)); if (present(idx)) CDF(:,idx) = C
                        PDF = 0.0; do k = 2,n; PDF(k) = (2.0/n)/(C(k+1)-C(k-1)); end do; if (oscale) PDF = PDF/maxval(PDF)
                        
                        ! this really should be replaced by a most likelihood fit of a gaussian distribution
                        avg = sum(C(2:n))/(n-1); var = sum((C(2:n)-avg)**2)/(n-2); X = (C-avg)/sqrt(2.0*var)
                        Y = (1.0 + erf(X))/2.0; Z = exp(-X**2); if (.not. oscale) Z = Z/sqrt(twopi*var)
                        
                        if (output$vis) then
#ifdef SILO
                                e = dbputcurve(db, STR("CDF_"//v), C, P, DB_DOUBLE, n+1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("PDF_"//v), C, PDF, DB_DOUBLE, n+1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("gaussian_CDF_"//v), C, Y, DB_DOUBLE, n+1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("gaussian_PDF_"//v), C, Z, DB_DOUBLE, n+1, DB_F77NULL, e)
#else
                                write (12,'(g)') "# CDF"
                                do k = 1,n+1; write (12,'(2g)') C(k), P(k); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PDF"
                                do k = 1,n+1; write (12,'(2g)') C(k), PDF(k); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# CDF [gaussian]"
                                do k = 1,n+1; write (12,'(2g)') C(k), Y(k); end do
                                write (12,'(g)') "", ""
                                
                                write (12,'(g)') "# PDF [gaussian]"
                                do k = 1,n+1; write (12,'(2g)') C(k), Z(k); end do
                                write (12,'(g)') "", ""
#endif
                        end if
                end if
                
#ifndef SILO
                if (output$vis) close (12)
#endif
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
                
                do k = 1,ns; write (31,'(32g)') t, (k-1)*dk, log10(PSD(k,1:idx)); end do
                write (31,'(g)') "", ""; flush(31)
        end if
        
        ! output distribution of values
        if (output$cdf) then
                inquire (32, opened=o)
                
                if (.not. o) then
                        open(32, file="CDF"); call head(32, (/"t", "percentile", DVAR(1:idx)/))
                end if
                
                do k = 1,n+1; write (32,'(32g)') t, real(k-1)/n, CDF(k,1:idx); end do
                write (32,'(g)') "", ""; flush(32)
        end if
end subroutine fflush

end
