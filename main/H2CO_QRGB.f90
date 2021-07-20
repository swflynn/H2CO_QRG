!=============================================================================80
!                       QRG grid and Collocation for H2CO                      !
!=============================================================================80
!       Discussion:
!Generate a quasiregular grid for H2CO PES and then use the collocation method
!to compute the spectra
!==============================================================================!
!       Modified:                                                              !
!   16 July 2021                                                               !
!       Authors:                                                               !
!   Shane Flynn and Vladimir Mandelshtam                                       !
!==============================================================================!
!CH2O of Carter 1997 Mol. Phy. The input distance is in bohr and angles are in
!radian, the return energy is in cm-1
!       O4
!       |
!       |
!       |
!       |1
!       C3
!      /  \
!     /    \
!    /      \
!  H1       H2
!the coordinates are polyspherical for bond vectors
!    r(1): C3-H1
!    r(2): C3-H2
!    r(3): C3-O4
!    r(4): th1
!    r(5): th2
!    r(6): phi
!the angles are defined based on the fact that C3-O4 is z axis,
!H1 is on x-z plane with positive projection
!==============================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!check unit conversion for grid point generation vs NC calc, make consistent
!store matricies to file
!store seed to file
!generate restart file
!make collocation a choice between pseudo/quasi random points
!FindRange has scaling in diag4, not in grid, need to decide if we keep, add flag
!Fix documentation for subroutines and main
!==============================================================================!
module QRGB_H2CO_grid_diag_mod
implicit none
!==============================================================================!
!NG              ==>Number of points to generate
!Dim_internal                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!E_cut                ==>Energy Cutoff Contour (1/cm input)
!==============================================================================!
double precision, parameter :: hbar = 1d0
double precision, parameter :: bohr = 0.52917721092
double precision, parameter :: autoeV = 27.211385
double precision, parameter :: eVtoau = 3.674932379D-2
double precision, parameter :: autocm = 2.194746313D5
double precision, parameter :: autokcalmol = 627.51
double precision, parameter :: Ktoau = 3.1668114D-6
double precision, parameter :: cmtoau = 4.5563352527D-6
double precision, parameter :: Hmass = 1837.15264758
double precision, parameter :: Cmass = 12.0000000 *Hmass/1.00782503223
double precision, parameter :: Omass = 15.99491461957 *Hmass/1.00782503223
double precision, parameter :: pi=acos(-1d0),d2r=pi/180d0
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
integer :: NG
integer(kind = 8) :: skip                                 !skip for sobol points
integer, parameter :: Natoms = 4, d = Natoms*3, d1 = d-6
double precision :: E_cut, Delta_E, rrange(d1), rrange1(d1), rrange2(d1)
double precision :: rmin(d1), rmax(d1), Pmax, gama
double precision, allocatable :: alpha(:), rg(:,:)
logical :: quasi
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!Compute mass of each atom, assumes water as input (H,O only)
!==============================================================================!
implicit none
character(len=2) :: atom
double precision :: Atom_Mass
!==============================================================================!
if (atom == 'H' .or. atom == 'h') then
  Atom_mass = Hmass
elseif (atom == 'O' .or. atom == 'o') then
  Atom_mass = Omass
elseif (atom == 'C' .or. atom == 'c') then
  Atom_mass = Cmass
else
  write(*,*) 'atom ', atom, ' is not recognized'
  STOP 'Check Function ==> Atom_Mass '
endif
end function Atom_Mass
!==============================================================================!
function P(r,Vr)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(r)
!r              ==>(d1) ith particles coordinate r^i_1,..,r^i_d
!gamma          ==>Scaling parameter for P: suggested d1/2
!==============================================================================!
implicit none
double precision :: r(d1), P, Vr
!==============================================================================!
call h2copot_carter97(r,Vr)                         ! H2CO potential units: 1/cm
!Vr=Vr*cmtoau     !diagonalization code uses this, also removed hardcoded cutoff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!come back to this i need to make it consistent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
P=1d-100
if(Vr<E_cut) P=((E_cut+Delta_E-Vr)*cmtoau)**(gama)
end function P
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!a              ==>uniform pseudo-random number
!==============================================================================!
implicit none
integer :: Nmin, Nmax, random_integer
double precision :: a
!==============================================================================!
call random_number(a)
random_integer = floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function Pair_rpl(r1,r2,sigma1,sigma2)
!==============================================================================!
!Pairwise repulsive energy between grid points
!==============================================================================!
!r1             ==>(d1) ith atoms coordinates
!r2             ==>(d1) jth atoms coordinates
!dist           ==>distance between points r1-r2
!==============================================================================!
implicit none
integer :: i
double precision :: r1(d1), r2(d1), Pair_rpl, sigma1, sigma2, dist
!==============================================================================!
dist = sqrt(sum(((r1(:)-r2(:))*rrange(:))**2))
Pair_rpl = (sigma1/dist)**(d1+9)+(sigma2/dist)**(d1+9)
end function Pair_rpl
!==============================================================================!
subroutine Metropolis_MC(N_MC,r1,Pr,Vr,step)
!==============================================================================!
!
!==============================================================================!
!r1             ==>Equilibrium geometry taken from Carter's 1997 paper
!               ==>
!               ==>
!==============================================================================!
implicit none
logical :: first
data first/.true./
integer :: N_MC, n, accept, my_size
double precision :: P_trial, Pr, a, Vr, V_trial, step
double precision :: r1(d1), Delr(d1), r_trial(d1)
integer, allocatable :: seed(:)
!==============================================================================!
if (first) then
  r1 = (/1.10064/bohr,1.10064/bohr,1.20296/bohr,121.65*d2r,121.65*d2r,180.*d2r/)
  Pr = P(r1,Vr)
  write(*,*) 'Initial geometry (internal coordinates) r0=', r1
  write(*,*) 'Energy(r0)=', Vr
!==============================================================================!
!initialize random_number seed
!Combe back will need to store the seed for restart calculations
!==============================================================================!
  call random_seed(size=my_size)
  allocate(seed(my_size))
  seed=my_size+1                         !gfortran seed must be larger than size
  call random_seed(put=seed)                      !seed must be an integer array
  first=.false.
  step=0.2
endif
accept=0
do n=1,N_MC
  call random_number(Delr)
  r_trial=r1+step*(2*Delr-1.0)
  P_trial=P(r_trial,V_trial)
  call random_number(a)
  if(P_trial/Pr > a) then
    r1=r_trial
    Pr=P_trial
    Vr=V_trial
    accept=accept+1
  endif
enddo
if(dble(accept)/N_MC > 0.6) step=step*1.1
if(dble(accept)/N_MC < 0.4) step=step*0.9
end subroutine Metropolis_MC
!==============================================================================!
subroutine histogram(r)
!==============================================================================!
!
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
implicit none
integer :: i,j,k
double precision, parameter :: rmin=0d0, rmar=2d0,dr=(rmar-rmin)/100
double precision :: hist(1:100), r(d1,NG), dist_aver, dist, dist1, Vr, Pr
!==============================================================================!
do i=1,NG
  Pr = P(r(:,i),Vr)
  dist1 = 1d10
  do j=1,NG
    if (i.ne.j) then
      dist = sqrt(sum((rrange(:)*(r(:,i)-r(:,j)))**2))
      if (dist < dist1) dist1 = dist
    endif
  enddo
  dist_aver=dist_aver+dist1*Pr**(1./d1)
enddo
dist_aver = dist_aver/NG
write(*,*) 'dist_aver =', dist_aver
hist=0d0
do i=1,NG
  Pr = P(r(:,i),Vr)
  do j=1,NG
    if (i.ne.j) then
      dist = sqrt(sum((rrange(:)*(r(:,i)-r(:,j)))**2))*Pr**(1./d1)/dist_aver
      k = (dist-rmin)/dr+1
      if (k>100) k = 100
      hist(k) = hist(k)+1
    endif
  enddo
enddo
open (8, file='hist.dat')
do i=1,100
  write(8,*) rmin+dr*(i-0.5),hist(i)
enddo
write(8,*)
end subroutine histogram
!==============================================================================!
subroutine find_range(r)
!==============================================================================!
!
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
implicit none
integer :: i, k
double precision :: rmin(d1), rmax(d1), r(d1,NG)
!==============================================================================!
rmin(:) = r(:,1)
rmax(:) = r(:,1)
do i=2,NG
  do k=1,d1
    if (rmin(k) > r(k,i)) rmin(k) = r(k,i)
    if (rmax(k) < r(k,i)) rmax(k) = r(k,i)
  enddo
enddo
do k=1,d1
  write(*,*) k,'  range=',rrange(k)
  rrange(k) = 1/(rmax(k)-rmin(k))
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! increase the range for the collocation points
!    rmin=rmin-0.2*rrange
!    rmax=rmax+0.2*rrange
!    rrange2=rmax-rmin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine find_range
!==============================================================================!
function cross(A,B)
!==============================================================================!
!
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
real (kind=8) :: A(3), B(3), cross(3)
!==============================================================================!
cross(1) = A(2)*B(3)-A(3)*B(2)
cross(2) = A(3)*B(1)-A(1)*B(3)
cross(3) = A(1)*B(2)-A(2)*B(1)
end function cross
!==============================================================================!
subroutine XYZms_to_Internal(x,r,flag)
!==============================================================================!
! mass-scaled XYZ to internal and back
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
logical :: flag
integer :: ix
real (kind=8) :: r(d1), r1, r2, r3, x(d), C(3), O(3), H1(3), H2(3), r1vec(3)
real (kind=8) :: r2vec(3), r3vec(3), r1r3normvec(3), r2r3normvec(3)
real (kind=8) :: costheta1, costheta2, cosphi, r1r3norm, r2r3norm, ksi
!==============================================================================!
if (flag) then                    !transform space fixed to internal coordinates
!==============================================================================!
  C(:) = x(1:3)/sqrt(Cmass)                 !unscale the mass-scaled coordinates
  O(:) = x(4:6)/sqrt(Omass)
  H1(:) = x(7:9)/sqrt(Hmass)
  H2(:) = x(10:12)/sqrt(Hmass)
  r1vec(:) = H1(:)-C(:)
  r1 = sqrt(sum(r1vec(:)**2))
  r2vec(:) = H2(:)-C(:)
  r2 = sqrt(sum(r2vec(:)**2))
  r3vec(:) = O(:)-C(:)
  r3 = sqrt(sum(r3vec(:)**2))
  costheta1 =dot_product(r1vec,r3vec)/(r1*r3)
  costheta2 =dot_product(r2vec,r3vec)/(r2*r3)
  r1r3normvec = cross(r1vec,r3vec)
  r1r3norm = sqrt(dot_product(r1r3normvec,r1r3normvec))
  r2r3normvec = cross(r2vec,r3vec)
  r2r3norm = sqrt(dot_product(r2r3normvec,r2r3normvec))
  cosphi = dot_product(r1r3normvec,r2r3normvec)/(r1r3norm*r2r3norm)
  r = (/r1,r2,r3,acos(costheta1),acos(costheta2),acos(cosphi)/)
  if(dot_product(r1r3normvec,r2vec)>0d0) r(6)=2*pi-r(6)
!==============================================================================!
else                        !transform internal to space fixed (XYZ) coordinates
!==============================================================================!
  x(:) = 0d0
  x(5) = r(3)                       !y(O)
  x(7) = -sin(r(4))*r(1)            !x(H1)
  x(8) = cos(r(4))*r(1)             !y(H1)
  x(11) = cos(r(5))*r(2)            !y(H2)
  ksi = sin(r(5))*r(2)              !distance from H2 to the y axis
  x(10) = ksi*cos(pi-r(6))          !x(H2)
  x(12) = ksi*sin(pi-r(6))          !z(H2)
  x(5) = x(5)*sqrt(Omass)           !mass-scale back
  x(7:12) = x(7:12)*sqrt(Hmass)     !mass-scale back
endif
!==============================================================================!
end subroutine XYZms_to_Internal
!==============================================================================!
function Phi(i,r)
!==============================================================================!
!d1-dimensional Gaussian
! Phi(i,r) = exp(alpha_i*(r_i-r)**2)
! r_ -- Gaussian centers (internal coordinates)
!               ==>
!               ==>
!==============================================================================!
integer :: i
real (kind=8) :: Phi, r(d1)
!==============================================================================!
Phi = exp(-alpha(i)*sum((rrange1(:)*(r(:)-rg(:,i)))**2))  !unnormalized Gaussian
end function Phi
!==============================================================================!
subroutine hpsort1(N,A,B,C)
!==============================================================================!
!Sort eigenvalues (real/complex)
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
Integer, Intent(In) :: N
Double Precision, Intent(Inout) :: A(N),B(N),C(N)
Integer :: I,IR,J,L
Double Precision :: RA,RB,RC
!==============================================================================!
L=N/2+1
IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
10  continue
    if (L > 1) then
       L=L-1
       RA=A(L)
       RB=B(L)
       RC=C(L)
    else
       RA=A(IR)
       RB=B(IR)
       RC=C(IR)
       A(IR)=A(1)
       B(IR)=B(1)
       C(IR)=C(1)
       IR=IR-1
       if (IR.eq.1) then
          A(1)=RA
          B(1)=RB
          C(1)=RC
          return
       endif
    endif
    I=L
    J=L+L
20  if (J.le.IR) then
       if (J < IR) then
          if (A(J) < A(J+1)) J=J+1
       endif
       if (RA < A(J)) then
          A(I)=A(J)
          B(I)=B(J)
          C(I)=C(J)
          I=J; J=J+J
       else
          J=IR+1
       endif
       goto 20
    endif
    A(I)=RA
    B(I)=RB
    C(I)=RC
    goto 10
end subroutine hpsort1
!==============================================================================!
subroutine generate_collocation_points(Nc,rc,V_)
!==============================================================================!
!generate a single collocation point r(d1), in the range
!pesudo is a global logical, if true use pseudo-random points, else quasi-random
!Uses the Sobol Sequence from sobol.f90 made available by John Burkardt
!https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html under GNU LGPL
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
use sobol       !for quasi-random sobol sequence
double precision :: rc(d1,Nc), dummy, V_(Nc), r(d1)
integer :: n, Nc
logical :: first
data first/.true./
!==============================================================================!
n=1
do while (n <= Nc)
!==============================================================================!
  if(quasi) then                 !Sobol Sequence from sobol.f90 by John Burkardt
    r = i8_sobol(int(d1, 8), skip)
  else
    call random_number(r)                                 !pseudo-random numbers
  endif
!==============================================================================!
  rc(:,n) = rmin + r * rrange2                             !Scale to span domain
  call random_number(dummy)
  if (P(rc(:,n),V_(n))/Pmax > dummy) then
    if (V_(n)<0d0) then
      write(*,*) 'n=', n, '  V=', V_(n)/cmtoau
      write(*,*) 'r=',rc(:,n)
      STOP 'Potential Energy < 0: Check Function ==> generate_collocation_points'
    endif
    if (first) write(3,*) V_(n)/cmtoau
    n=n+1
  end if
end do
if(first) close(3)
first=.false.
end subroutine generate_collocation_points
!==============================================================================!
subroutine matrix_elements(S,H,flag)
!===========================================================================!
!flag=.true. ---> Using two grids   flag=.false. ---> single grid rc=rg
!Construct S and H using the Gaussian centers rg(d1,NG) and collocation points rc(d1,NG)
!Using 7-point stencil (finite difference) to compute the laplassian of Phi
!==============================================================================!
!               ==>
!               ==>
!==============================================================================!
integer :: i, n, k, l
double precision :: S(NG,NG), H(NG,NG), V_(NG), rc(d1,NG)
logical :: flag
real (kind=8) :: y(d),yk,rr(d1)
real (kind=8) :: HH
real (kind=8), parameter :: Delta_Y = 0.1
integer, parameter :: stencil=1
real(kind=8) :: a(-stencil:stencil) = (/-1.,d*2.,-1./)/(2*delta_Y**2)
!==============================================================================!
if (flag) then          !additional collocation points and compute the potential
  call generate_collocation_points(NG,rc,V_)
else                           !Use the Gaussian ceneters for collocation points
  do n=1,NG
    rc(:,n) = rg(:,n)
    call h2copot_carter97(rc(:,n),V_(n))             !H2CO potential units: 1/cm
    V_(n) = V_(n)*cmtoau
    write(3,*) V_(n)/cmtoau
  enddo
  write(3,*)
  flush(3)
endif
S=0d0
H=0d0
do i=1,NG
  do n=1,NG
    rr = rc(:,n)
    S(n,i) = Phi(i,rr)
    HH = (V_(n)+a(0))*S(n,i)
    call XYZms_to_Internal(y,rr,.false.)                                 !r--->y
    do k=1,d
      yk = y(k)
      do l=-stencil,stencil
        if (l.ne.0) then
          y(k) = yk+l*Delta_Y
          call XYZms_to_Internal(y,rr,.true.)                  !y-->r (unscaled)
          HH = HH+a(l)*Phi(i,rr)
        endif
      enddo
      y(k) = yk
    enddo
    H(n,i) = HH
  enddo
enddo
end subroutine matrix_elements
!==============================================================================!
end module QRGB_H2CO_grid_diag_mod
!==============================================================================!
program main
use QRGB_H2CO_grid_diag_mod
implicit none
!==============================================================================!
!               ==>
!               ==>
!               ==>
!               ==>
!               ==>
!               ==>
!               ==>
!alpha0         ==>Parameter to scale Gaussian Widths
!S              ==>(NG,NG) Overlap Matrix
!H              ==>(NG,NG) Hamiltonian Matrix
!==============================================================================!
integer :: i, j, k, n, N_MC_QRG, accept, Ndiag, Nblocks
double precision :: Vpot, time0, time1, time2, delE, dummy, mv_cutoff, P_trial
double precision :: step, Vr, sigma_trial, E_total, V_trial, Pr, dist, dist1
double precision :: dist_aver,step,Vr, hist(1:100), delr(d1), r_trial(d1)
double precision :: rmin(d1), rmax(d1), x2, alpha0, V0
double precision, allocatable :: r0(:), r(:,:), sigma(:), V_(:), S(:,:), H(:,:)
double precision, allocatable :: S_save(:,:), H_save(:,:)
character(len=50) :: Gaussian_centers, Collocation_points
!=============================================================================!
!                    LLAPACK Diagonalization Variables
!=============================================================================!
integer :: info,Lwork
double precision, allocatable :: work(:), Er(:), Ei(:), BETA(:), VL(:,:), VR(:,:)
!=============================================================================!
!                             Read Input Data File
!=============================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!==============================================================================!
program main_grid
  use  QRGB_H2CO_grid_mod
  !==============================================================================!

  !==============================================================================!
  !                              Read Input File                                 !
  !==============================================================================!
  call cpu_time(time1)
  read(*,*) NG
  read(*,*) E_cut, Delta_E
  read(*,*) N_MC_QRG

   !==============================================================================!
  allocate(r(d1,NG),sigma(NG),V_(NG))

  !==============================================================================!
  !                   Generate initial random grid
  !==============================================================================!
  open(123,file='Metropolis.dat')
  write(123,*) '# Equilibration'
  do i=1,100
     call Metropolis_MC(100,r(:,1),Pr,Vr,step)
     write(123,*) i,Vr, step
  enddo


  write(123,*) '# Start the grid'
  do i=2,NG
     r(:,i)=r(:,i-1)
     call Metropolis_MC(1000,r(:,i),Pr,Vr,step)
     write(123,*) i, Vr, step
  enddo
  close(123)

  open(21,file='grid_rand.dat')
  do i=1,NG
     write(21,*) r(:,i)
  enddo
  close(21)
  !==============================================================================!
  !                         Generate QRG (greedy search)
  !==============================================================================!
  open(72,File='Simul_anneal.dat')
  call find_range(r)
  do i=1,NG
     sigma(i)=P(r(:,i),V_(i))**(-1./d1)
  enddo
  open(124,file='potential.dat')
  write(124,*) '# random grid'
  do i=1,NG
     write(124,*) i,V_(i)
  enddo
  write(124,*)
  write(124,*) '# QRG grid'
  E_total=0
  do i=2,NG
     do j=1,i-1
        E_total=E_total+Pair_rpl(r(:,i),r(:,j),sigma(i),sigma(j))
     enddo
  enddo

  write(*,*) 'Generated initial random grid, E_total = ',E_total


  write(*,*) 'Start optimizing QRG'

  accept=0
  step=0.01
  do n=1,NG*N_MC_QRG
     k=random_integer(1,NG)                               !Select Atom to Move
     call random_number(delr)
     r_trial=r(:,k)+step*(2*delr-1)
     P_trial=P(r_trial,V_trial)
     if(P_trial > 1d-20) then                         !Only consider if V(trial)<Ecut
        sigma_trial=P_trial**(-1./d1)
        delE=0d0
        do j=1,NG
           if(j.ne.k) delE=delE+Pair_rpl(r(:,j),r_trial,sigma(j),sigma_trial) &
                   -Pair_rpl(r(:,j),r(:,k),sigma(j),sigma(k))     !Energy change trial move
        enddo
        if(delE <= 0d0) then
           r(:,k)=r_trial(:)
           sigma(k)=sigma_trial
           V_(k)=V_trial
           accept=accept+1
           E_total=E_total+delE
        endif
     endif
     !for MMC want acceptance ~40-60%, adjust trial movement displacement accordingly
     if(mod(n,NG*10)==0)then
        if(accept/NG*10. < 0.4) step=step*0.9
        if(accept/NG*10. > 0.6) step=step*1.1
        accept=0
        E_total=0
        do i=2,NG
           do j=1,i-1
              E_total=E_total+Pair_rpl(r(:,i),r(:,j),sigma(i),sigma(j))
           enddo
        enddo
        write(72,*) n,E_total,step
     endif
  enddo
  close(72)
  write(*,*) 'Etotal after minimization=', E_total
  open(22,File='grid_QRG.dat')
  do i=1,NG
     write(124,*) i,V_(i)
     write(22,*) r(:,i)
  enddo
  close(124)
  close(22)
  write(*,*) 'Compute the pair correlation function'
  open (8, file='hist.dat')
  call histogram(r)
  close(8)
  call cpu_time(time2)
  write(*,*) 'Simulation Time==> ',time2-time1

end program main_grid

!=============================================================================!
!=============================================================================!
!=============================================================================!
!=============================================================================!
!End of grid code
!start of diag code
!=============================================================================!
!=============================================================================!
!=============================================================================!
!=============================================================================!

program QRGB_H2CO_diag

  !=============================================================================!
  !                             Read Input Data File
  !=============================================================================!
  call cpu_time(time0)
  time1=time0

  read(*,*) Gaussian_centers, NG
    read(*,*) Ndiag, Nblocks
  ! Nblocks --- Nblocks*NG = number of collocation points added for each diagonalization
  ! Ndiag ---  number of diagonalizations
  read(*,*) alpha0
  read(*,*) E_cut
  read(*,*) Delta_E
  E_cut=E_cut*cmtoau
  Delta_E=Delta_E*cmtoau
  read(*,*) quasi
!  read(*,*) Delta_Y
  skip=NG   !set random number generator skip
  !=============================================================================!
  !                               Allocations
  !=============================================================================!
  allocate(rg(d1,NG),alpha(NG),S(NG,NG),H(NG,NG),r0(d1))
  !=============================================================================!
  !                           Read GridPoints x(d,NG)
  !=============================================================================!
  write(*,*) 'Number of Gaussains,  NG =', NG
  write(*,*) ' alpha0=',alpha0
  write(*,*) 'Number of collocation points added for each diagonalization = Nblocks * NG', Nblocks,'*',NG
  write(*,*) 'Total number of diagonalizations =', Ndiag
  write(*,*) 'File with Gaussian centers =', Gaussian_centers
  open(1,file=Gaussian_centers)
  do i=1,NG
     read(1,*) rg(:,i)
  enddo
  close(1)
  call find_range(NG,rg)
  do i=1,NG
     alpha(i)=1d20                            !large distance for placeholder
     do j=1,NG
        if(j.ne.i) then
!           x2=sum((rg(:,i)-rg(:,j))**2)               !distance between gridpoints
           x2=sum((rrange1(:)*(rg(:,i)-rg(:,j)))**2)    !scaled distance between gridpoints
           if(x2<alpha(i)) alpha(i)=x2
        endif
     enddo
     alpha(i)=alpha0/alpha(i)
  enddo
  open(unit=21,file='eig.dat')
  open(unit=22,file='freq.dat')
  open(3,file='potential.dat')
  write(21,*) '# NG=',NG,' alpha0=',alpha0
  write(22,*) '# NG=',NG,' alpha0=',alpha0
  r0 =(/1.10064/bohr,1.10064/bohr,1.20296/bohr,121.65*d2r,121.65*d2r,180.*d2r/)
!  r0 =(/2.03832, 2.03832, 2.43021, 119.697, 119.697, 180.000/)
  Pmax=P(r0,V0)
  write(*,*) 'Pmax ==> ', Pmax
  write(*,*) 'V0 ==> ', V0
  write(*,*) 'r0 ==> ', r0

  write(*,*) 'First diagonalization Round. Using Gaussian centers for collocation points'
  call matrix_elements(S,H,.false.)
  if(Ndiag > 0) then
     allocate(S_save(NG,NG),H_save(NG,NG))
     S_save=S
     H_save=H
  endif
  do k=0,Ndiag
     write(*,*) 'Diagonalization Round =', k, ' Totalal number of collocation points =', k*Nblocks*NG
     if(k > 0) then
        do n=1,Nblocks
           if(k==1.and.n==1) then
              S=S_save
              H=H_save
           else
              call matrix_elements(S,H,.true.)
           endif
           do i=1,NG
              do j=1,NG
                 H_save(j,i)=H_save(j,i)+dot_product(H(:,i),S(:,j))
                 S_save(j,i)=S_save(j,i)+dot_product(S(:,i),S(:,j))
              enddo
           enddo
        enddo
     endif
     call cpu_time(time2)
     write(*,*) 'Constructed the matrices time ==>', time2-time1
     time1=time2
     write(*,*)  "Solving the generalized eigenvalue problem..."
     flush(6)

     Lwork=max(1,8*NG)
     allocate(work(Lwork),Er(NG),Ei(NG),BETA(NG))
     if(k > 0) then
        S=S_save
        H=H_save
     endif
     call dggev('n','n',NG,H,NG,S,NG,Er,Ei,BETA,VL,NG,VR,NG,work,Lwork,info)
     Er(:)=Er(:)/BETA(:)
     Ei(:)=Ei(:)/BETA(:)
     call hpsort1(NG,Er,Ei,BETA)
     write(*,*) 'info ==> ', info
     do i=1,50
        write(21,*) k*Nblocks*NG, Er(i)/cmtoau, Ei(i)/cmtoau, BETA(i)/cmtoau
        write(22,*) k*Nblocks*NG, (Er(i)-Er(1))/cmtoau
     enddo
     deallocate(work,Er,Ei,BETA)
     write(21,*)
     write(22,*)
     flush(21)
     flush(22)

     call cpu_time(time2)
     write(*,*) 'Diagonalization time ==>', time2-time1
     time1=time2
     flush(6)

  enddo
  write(*,*) 'Total time ==>', time2-time0

end program QRGB_H2CO_diag

