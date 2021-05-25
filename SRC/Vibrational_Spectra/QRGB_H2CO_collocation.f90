!==============================================================================!
!              H2CO Vibrational Energy Calculations (Collocation)              !
!==============================================================================!
!This implementation is adapted from Manzhos and Carrington 2016 paper
!"Using an internal coordinate Gaussian basis and a space-fixed Cartesian
!coordinate kinetic energy operator to compute a vibrational spectrum with
!rectangular collocation"
!Authors matlab code is available in their SI and was used to develop this code
!==============================================================================!
!       Modified:
!   19 May 2021
!       Authors:
!   Shane Flynn and Vladimir Mandelshtam
!==============================================================================!
! TO DO: finish documentation for subroutines etc
!==============================================================================!
!The following PES is used for this project:              (h2copot_carter97.f90)
!==============================================================================!
!CH2O of Carter 1997 Mol. Phy.
!The coordinates are polyspherical for bond vectors
!
!The input distance is in bohr and angles are in radian, and the return energy
!is in cm-1
! Carter's original pes : hartree/rad^n
!       O4
!       |
!       |
!       |
!       |1
!       C3
!      /  \
!     /     \
!    /        \
!  H1         H2
!r(1): C3-H1
!r(2): C3-H2
!r(3): C3-O4
!r(4): th1
!r(5): th2
!r(6): phi
!the angles are defined based on the fact that C3-O4 is z axis,
!H1 is on x-z plane with positive projection
!==============================================================================!
module QRGB_H2CO_collocation_mod
implicit none
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
double precision, parameter :: Dmass = 2.01410177812 * Hmass/1.00782503223
double precision, parameter :: Cmass = 12.0000000 * Hmass/1.00782503223
double precision, parameter :: Omass = 15.99491461957 * Hmass/1.00782503223
double precision, parameter :: Fmass = 18.99840316273 * Hmass/1.00782503223
double precision, parameter :: Clmass = 34.968852682 * Hmass/1.00782503223
double precision, parameter :: Brmass = 78.9183376 * Hmass/1.00782503223
double precision, parameter :: Imass = 126.9044719 * Hmass/1.00782503223
double precision, parameter :: pi = acos(-1d0), d2r = pi/180d0
!=============================================================================!
!                            Global Variables                                 !
!=============================================================================!
integer, parameter :: Natoms = 4, d = Natoms*3, d1 = d-6
integer :: NG, NC, Nthreads
double precision,allocatable :: alpha(:), x1(:,:), r(:,:), V_(:)
!==============================================================================!
!Natoms             ==>Number of atoms in the system
!d                  ==>Full system dimensionality (Natoms*(x+y+z))
!d1                 ==>System subspace (depends on the molecule)
!NG                 ==>Number of Gaussian Basis Functions
!NC                 ==>Number of Collocation Points
!alpha()            ==>Gaussian Widths (small alpha=broad Gaussian)
!x1(d1,NC)          ==>Gaussian Centers
!r(d1,NC)           ==>Subspace coordinates
!V_(NC)             ==>Potential at each NC
!==============================================================================!
contains
!==============================================================================!
function cross(A, B)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
implicit none
double precision :: A(3), B(3), cross(3)
!==============================================================================!
cross(1)=A(2)*B(3)-A(3)*B(2)
cross(2)=A(3)*B(1)-A(1)*B(3)
cross(3)=A(1)*B(2)-A(2)*B(1)
end function cross
!==============================================================================!
subroutine XYZms_to_Internal(x, r, flag)
!==============================================================================!
!mass-scaled XYZ to internal coordinates and back
!==============================================================================!
!
!
!==============================================================================!
implicit none
logical :: flag
double precision :: r(d1), x(d)
double precision :: C(3), O(3), H1(3), H2(3), r1vec(3), r2vec(3), r3vec(3)
double precision :: r1r3normvec(3), r2r3normvec(3), r1, r2, r3, costheta1
double precision :: costheta2, cosphi, r1r3norm, r2r3norm, ksi
!==============================================================================!
if(flag) then !transform space fixed to internal coordinates (Taken from Tucker)
  C(:) = x(1:3)/sqrt(Cmass)                     !unscale mass-scaled coordinates
  O(:) = x(4:6)/sqrt(Omass)
  H1(:) = x(7:9)/sqrt(Hmass)
  H2(:) = x(10:12)/sqrt(Hmass)
!==============================================================================!
  r1vec(:) = H1(:)-C(:)
  r1 = sqrt(sum(r1vec(:)**2))
  r2vec(:) = H2(:)-C(:)
  r2 = sqrt(sum(r2vec(:)**2))
  r3vec(:) = O(:)-C(:)
  r3 = sqrt(sum(r3vec(:)**2))
!==============================================================================!
  costheta1 = dot_product(r1vec, r3vec)/(r1*r3)
  costheta2 = dot_product(r2vec, r3vec)/(r2*r3)
  r1r3normvec = cross(r1vec, r3vec)
  r1r3norm = sqrt(dot_product(r1r3normvec, r1r3normvec))
  r2r3normvec = cross(r2vec, r3vec)
  r2r3norm = sqrt(dot_product(r2r3normvec, r2r3normvec))
  cosphi = dot_product(r1r3normvec, r2r3normvec)/(r1r3norm*r2r3norm)
!==============================================================================!
  r = (/r1, r2, r3, acos(costheta1), acos(costheta2), acos(cosphi)/)    !correct
!==============================================================================!
  if(dot_product(r1r3normvec, r2vec) > 0d0) r(6) = 2*pi-r(6)
!==============================================================================!
else                        !transform internal to space fixed (XYZ) coordinates
  x(:) = 0d0
  x(5) = r(3)                     !y(O)
  x(7) = -sin(r(4))*r(1)          !x(H1)
  x(8) = cos(r(4))*r(1)           !y(H1)
  x(11) = cos(r(5))*r(2)          !y(H2)
  ksi = sin(r(5))*r(2)            !distance from H2 to the y axis
  x(10) = ksi*cos(pi-r(6))        !x(H2)
  x(12) = ksi*sin(pi-r(6))        !z(H2)            !!!!transformation ends here
!==============================================================================!
!!!!  transformation ends here
  x(5) = x(5)*sqrt(Omass)
  x(7:12) = x(7:12)*sqrt(Hmass)                                 !mass-scale back
endif
end subroutine XYZms_to_Internal
!==============================================================================!
function Phi(i, x, r1, flag)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
! d1-dimensional Gaussian
! Phi(i,x) = exp(alpha_i*(r_i-r1(x))**2)
! r_i -- Gaussian centers
! x(1:d) - collocation point (using cartesian mass-scaled coordinates)
! r1(1:d1) = r1(x) - internal coordinates
! flag = .true. ---> transfrom x-->r
implicit none
logical :: flag
integer :: i
double precision :: Phi, x(d), r1(d1)
!==============================================================================!
if(flag) call XYZms_to_Internal(x, r1, .true.)     !transform to internal coords
Phi = exp(-alpha(i)*sum((r1(:)-r(:,i))**2))               !unnormalized Gaussian
end function Phi
!==============================================================================!
subroutine hpsort(N,A,B)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
implicit None
integer, intent(in) :: N
double precision, intent(inout) :: A(N), B(N)
integer :: I, IR, J, L
double precision :: RA, RB
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
else
  RA=A(IR)
  RB=B(IR)
  A(IR)=A(1)
  B(IR)=B(1)
  IR=IR-1
  if (IR.eq.1) then
    A(1)=RA
    B(1)=RB
    return
  endif
endif
I=L
J=L+L
20 if (J.le.IR) then
    if (J < IR) then
      if (A(J) < A(J+1)) J=J+1
    endif
    if (RA < A(J)) then
      A(I)=A(J)
      B(I)=B(J)
      I=J; J=J+J
    else
      J=IR+1
    endif
  goto 20
endif
A(I)=RA
B(I)=RB
goto 10
end subroutine hpsort
!==============================================================================!
subroutine hpsort1(N,A,B,C)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
implicit None
integer, intent(in) :: N
double precision, intent(inout) :: A(N),B(N),C(N)
integer :: I,IR,J,L
double precision :: RA,RB,RC
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
function V(x, rr, flag)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
! This functions computes H2CO potential (atomic units)
! If flag=.true. it takes mass-scaled Cartesian coordinates x(1:12) of H2CO
! and transfers them into the 6 internal coordinates rr(1:6)
implicit none
logical :: flag
double precision :: Vr, V, rr(6), x(12), C(3), O(3), H1(3), H2(3), r1vec(3)
double precision :: r2vec(3), r3vec(3), r1r3normvec(3), r2r3normvec(3), r1, r2
double precision :: r3, costheta1, costheta2, cosphi, r1r3norm, r2r3norm, ksi
!==============================================================================!
if(flag) then   !transform x(1:d) --> rr(1:d1)    (transformation from Tucker's)
  C(:) = x(1:3)/sqrt(Cmass)                     !unscale mass-scaled coordinates
  O(:) = x(4:6)/sqrt(Omass)
  H1(:) = x(7:9)/sqrt(Hmass)
  H2(:) = x(10:12)/sqrt(Hmass)
!==============================================================================!
  r1vec(:) = H1(:)-C(:)
  r1 = sqrt(sum(r1vec(:)**2))
  r2vec(:) = H2(:)-C(:)
  r2 = sqrt(sum(r2vec(:)**2))
  r3vec(:) = O(:)-C(:)
  r3 = sqrt(sum(r3vec(:)**2))
!==============================================================================!
  costheta1 = dot_product(r1vec, r3vec)/(r1*r3)
  costheta2 = dot_product(r2vec, r3vec)/(r2*r3)
  r1r3normvec = cross(r1vec, r3vec)     !sine of cross doesnt make a diff. for V
  r1r3norm = sqrt(dot_product(r1r3normvec, r1r3normvec))
  r2r3normvec = cross(r2vec, r3vec)
  r2r3norm = sqrt(dot_product(r2r3normvec, r2r3normvec))
  cosphi = dot_product(r1r3normvec, r2r3normvec)/(r1r3norm*r2r3norm)
!==============================================================================!
  rr = (/r1, r2, r3, acos(costheta1), acos(costheta2), acos(cosphi)/)   !correct
!==============================================================================!
  if(dot_product(r1r3normvec, r2vec)>0d0) rr(6)=2*pi-rr(6)   !End transformation
endif
!==============================================================================!
call h2copot_carter97(rr, Vr)                                    !H2CO potential
V=Vr*cmtoau                                                 !now V is in Hartree
end function V
!==============================================================================!
subroutine matrix_elements_col_sq(Smat, Hmat)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
implicit none
integer :: i, j, n, k, l
double precision :: Smat(NG,NC),Hmat(NG,NC),y(d),yk,rr(d1)
!Using 5-point stencil (finite difference) to compute the laplassian of Phi
double precision, parameter :: ss = 1d-3                          !need to check
double precision, parameter :: a(-2:2) = (/1.,-16.,30.,-16.,1./)/(24*ss**2)
!==============================================================================!
do i = 1,NG
  do n = 1,Nc
    y(:)=(/0d0, 0d0, 0d0, 0d0, x1(1,n), 0d0, x1(2,n), x1(3,n), 0d0, x1(4,n), &
    x1(5,n), x1(6,n)/)
    Smat(i,n) = Phi(i, y, r(:,n), .false.)
    Hmat(i,n) = (V_(n)+d*a(0))*Smat(i,n)
    do k = 1,d
      do l = -2,2
        if(l.ne.0) then
          yk = y(k)
          y(k) = y(k)+l*ss
          Hmat(i,n) = Hmat(i,n)+a(l)*Phi(i, y, rr, .true.)
          y(k) = yk
        endif
      enddo
    enddo
  enddo
enddo
end subroutine matrix_elements_col_sq
!==============================================================================!
subroutine matrix_elements_col_rec(Smat, Hmat)
!==============================================================================!
!
!==============================================================================!
!
!
!==============================================================================!
implicit none
integer :: i, j, n, k, l, m, n1
double precision :: Smat(NG,NG), Hmat(NG,NG), y(d), yk, rr(d1), S(NG,NG), H_i(NG)
!Using 5-point stencil (finite difference) to compute the laplassian of Phi
double precision, parameter :: ss = 1d-3                          !need to check
double precision, parameter :: a(-2:2) = (/1.,-16.,30.,-16.,1./)/(24*ss**2)
!==============================================================================!
do m = 0,Nc/NG-1
  !$OMP parallel default(shared) private(i,j,n,n1,y,k,l,yk,rr)
  !$OMP do
  do i = 1,NG
    do n = 1,NG
      n1 = m*NG+n
      S(n,i) = Phi(i, y, r(:,n1), .false.)
    enddo
  enddo
  !$OMP enddo
!==============================================================================!
  do i=1,NG
    !$OMP do
    do n=1,NG
      n1=m*NG+n
      H_i(n) = (V_(n1)+d*a(0))*S(n,i)
      y(:)=(/0d0, 0d0, 0d0, 0d0, x1(1,n1), 0d0, x1(2,n1), x1(3,n1), 0d0, &
      x1(4,n1), x1(5,n1), x1(6,n1)/)
      do k = 1,d
        do l = -2,2
          if(l.ne.0) then
            yk = y(k)
            y(k) = y(k)+l*ss
            H_i(n) = H_i(n)+a(l)*Phi(i, y, rr, .true.)
            y(k) = yk
          endif
        enddo
      enddo
    enddo
    !$OMP enddo
!==============================================================================!
!                 Update the square Hmat and Smat matrices
!==============================================================================!
    !$OMP do
    do j = 1,NG
      Hmat(j,i) = Hmat(j,i)+dot_product(H_i, S(:,j))
      Smat(j,i) = Smat(j,i)+dot_product(S(:,i), S(:,j))
    enddo
    !$OMP enddo
  enddo
  !$OMP end parallel
enddo
end subroutine matrix_elements_col_rec
!==============================================================================!
end module QRGB_H2CO_collocation_mod
!==============================================================================!
program QRGB_H2CO_collocation
!==============================================================================!
use QRGB_H2CO_collocation_mod
!==============================================================================!
!NG             ==>Number of Gaussian Basis Functions (gridpoints)
!alpha0         ==>Flat Scaling Parameter for Gaussian Widths
!Smat           ==>(NG,NG) Overlap Matrix
!Hmat           ==>(NG,NG) Hamiltonian Matrix
!==============================================================================!
implicit none
character(len=50) :: Gaussian_centers, Collocation_points
integer :: i, j ,k, Nalpha
double precision :: x2, alpha0, alpha1, alphamin, alphamax
double precision, allocatable :: Smat(:,:), Hmat(:,:)
!==============================================================================!
!                            LLAPACK  variables
!==============================================================================!
integer :: info, Lwork
double precision, allocatable :: work(:), Er(:), Ei(:), BETA(:), VL(:,:)
double precision, allocatable :: VR(:,:), y(:)
!==============================================================================!
!                             Read Input Data File
!==============================================================================!
read(*,*) Gaussian_centers, NG
read(*,*) Collocation_points, Nc
if(mod(Nc,NG).ne.0) STOP 'Nc should be multiple of NG'
Nc = Nc+NG
if(Collocation_points == 'none') Nc=NG
read(*,*) alphamin, alphamax, Nalpha
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x1(d1,NC), r(d1,NC), V_(NC), alpha(NG), Smat(NG,NG), Hmat(NG,NG), y(d))
!==============================================================================!
!                           Read GridPoints x(d,NG)
!==============================================================================!
open(1, File = Gaussian_centers)
do i = 1,NG
  read(1,*) x1(:,i)
enddo
close(1)
write(*,*) 'Using collocation method  NG NC=',NG, NC
write(*,*) 'File with Gaussian centers=', Gaussian_centers
if(NC>NG) then  !add collocation points
  open(1, File = Collocation_points)
  write(*,*) 'File with additional collocation points== ', Collocation_points
  do i = NG+1,Nc
    read(1,*) x1(:,i)
  enddo
  close(1)
endif
open(3, File = 'potential.dat')
do i = 1,Nc
  y(:) = (/0d0, 0d0, 0d0, 0d0, x1(1,i), 0d0, x1(2,i), x1(3,i), 0d0, x1(4,i), &
  x1(5,i), x1(6,i)/)
  V_(i) = V(y, r(:,i), .true.)                    !H2CO potential (atomic units)
  write(3,*) V_(i)/cmtoau
enddo
close(3)
!==============================================================================!
!         Symmetric gaussians. Use nearest neighbor to determine alpha
!==============================================================================!
do i = 1,NG
  alpha(i) = 1d20                            !large distance for placeholder
  do j = 1,NG
    if(j.ne.i) then
      x2 = sum((r(:,i)-r(:,j))**2)          !distance between gridpoints
      if(x2<alpha(i)) alpha(i) = x2
    endif
  enddo
  alpha(i) = alphamin/alpha(i)
enddo
!==============================================================================!
alpha1 = alphamin
open(unit=21,File = 'alpha_eig.dat')
open(unit=22,File = 'alpha_freq.dat')
do k = 0,Nalpha
  alpha0 = alphamin+dble(k)*(alphamax-alphamin)/Nalpha
  alpha = alpha*alpha0/alpha1
  alpha1 = alpha0
  write(*,*) 'alpha0= ', alpha0 , ' Constructing Hmat and Smat...'
  if(Nc==NG) then
    call  matrix_elements_col_sq(Smat, Hmat)
  else
    call  matrix_elements_col_rec(Smat, Hmat)
  endif
  write(*,*)  "Solving the generalized eigenvalue problem..."
  flush(6)
!==============================================================================!
  Lwork = max(1,8*NG)
  allocate(work(Lwork), Er(NG), Ei(NG), BETA(NG))
  call dggev('n', 'n', NG, Hmat, NG, Smat, NG, Er, Ei, BETA, VL, NG, VR, NG, &
  work,Lwork,info)
  Er(:) = Er(:)/BETA(:)
  Ei(:) = Ei(:)/BETA(:)
  call hpsort1(NG, Er, Ei, BETA)
  write(*,*) 'info ==> ', info
!==============================================================================!
  do i = 1,50
    write(21,*) alpha0, Er(i)/cmtoau, Ei(i)/cmtoau, BETA(i)
    write(22,*) alpha0, Er(i)/cmtoau -  Er(1)/cmtoau
  enddo
  deallocate(work, Er, Ei, BETA)
  write(21,*)
  write(22,*)
  flush(21)
  flush(22)
!==============================================================================!
  flush(6)
enddo
end program QRGB_H2CO_collocation
