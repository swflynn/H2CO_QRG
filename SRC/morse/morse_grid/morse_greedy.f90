!=============================================================================80
!                       QRG Morse Implementation                  !
!==============================================================================!
!simple ndmorse implementation for testing greedy minimization
!implement move cutoff that scales wrt sigma.
!This is the most recent version with different integral_P options
!==============================================================================!
!    Modified:
!4 November 2020
!    Author:
!Shane Flynn
!==============================================================================!
module morse_mod
implicit none
!==============================================================================!
!Npoints              ==>Number of points to generate
!d                    ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!r(d2,Npoints)        ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!Uij(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                  ==>Potential Energy evaluation V(x_i)
!E_cut                ==>Distribution cutoff contour (Kcal/mol input)
!rmin(d)              ==>Minimum of normalization box size
!rmax(d)              ==>Maximum of normalization box size
!N_MMC_box            ==>Number of MMC Iterations to determine box-size
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid           ==>Number of MMC Iterations to optimize QRG
!MMC_freq             ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P           ==>Normalization constant for the distribtion P(x)
!c_LJ                 ==>Parameter for q-LJ pseudo-potential
!E_cut                ==>Energy Cutoff Contour (kcal/mol input)
!omega                ==>(d) Parameter for Morse Potential
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
integer::d,Npoints
double precision,allocatable,dimension(:)::omega
double precision,allocatable,dimension(:,:)::x
double precision::E_cut,integral_P,c_LJ
double precision::gammas
!==============================================================================!
contains
!==============================================================================!
function P_i(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none
double precision::x(d),P_i
!if(V(x)<E_cut) P_i=(E_cut-V(x))**(d/2.)/integral_P
if(V(x)<E_cut) P_i=(E_cut-V(x))**(gammas)/integral_P
if(V(x)>=E_cut) P_i=1d-20
end function P_i
!==============================================================================!
function V(x)
!==============================================================================!
!Potential Energy (sum of 1D morse potentials)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
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
integer::Nmin,Nmax,random_integer
double precision::a
!==============================================================================!
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
function Pair_rpl(x1,x2,sigma1,sigma2)
!==============================================================================!
!Pairwise repulsive energy between grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!==============================================================================!
implicit none
double precision::x1(d),x2(d),Pair_rpl,sigma1,sigma2,dist
dist=sqrt(sum((x1(:)-x2(:))**2))
!sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)
!sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)
Pair_rpl=(sigma1/dist)**(d+9)+(sigma2/dist)**(d+9)
end function Pair_rpl
!==============================================================================!
function Pair_LJ_NRG(x1,x2,sigma1,sigma2)
!==============================================================================!
!quasi-Lennard Jones pairwise energy between grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate q-LJ
!sigma          ==>Gaussian widths (c*sigma(P))
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2
a=sqrt(sum((x1(:)-x2(:))**2))
!sigma1=c_LJ*(P(x1)*Npoints)**(-1./d)
!sigma2=c_LJ*(P(x2)*Npoints)**(-1./d)
b=sigma2/a
a=sigma1/a
Pair_LJ_NRG=a**(d+9)-a**(d+3)+b**(d+9)-b**(d+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine domain_size_P(xmin,xmax,N_MMC_box)
!==============================================================================!
!Determine the box size for normalizing P(r) using Metropolis Monte Carlo
!This subroutine Computes rmin(d) and rmax(d)
!Be careful about making mv_cutoff large with MMC, PES does not go to infinity
!the surface may not have a well-defined minimum in higher dimensions
!==============================================================================!
!xmin(d)              ==>Minimum of normalization domain
!xmax(d)              ==>Maximum of normalization domain
!N_MMC_box            ==>Number of MMC Iterations to determine box-size
!==============================================================================!
double precision,parameter::mv_cutoff=0.1
integer::N_MMC_box,my_size,i,j
integer,allocatable::seed(:)
double precision::x_trial(d),x_i(d),s(d),xmin(d),xmax(d),V_trial
!==============================================================================!
!Set seed for random number generator (make it reproducible)
!==============================================================================!
call random_seed(size=my_size)
allocate(seed(my_size))
seed=my_size+1                           !gfortran seed must be larger than size
call random_seed(put=seed)                        !seed must be an integer array
!==============================================================================!
integral_P=1d0                             !Initially set to 1 so you can call P
x_i=0d0
xmin=x_i
xmax=x_i
do i=1,N_MMC_box
  call random_number(s)
  x_trial=x_i+mv_cutoff*(2*s-1)                !make trial move (-1,1) not (0,1)
  V_trial=V(x_trial)
  if(V_trial<E_cut) then
     x_i=x_trial
     do j=1,d
        if(xmin(j).gt.x_i(j)) xmin(j)=x_i(j)
        if(xmax(j).lt.x_i(j)) xmax(j)=x_i(j)
     enddo
  endif
enddo
end subroutine domain_size_P
!==============================================================================!
subroutine normalize_P(norm_method,xmin,xmax,Nevals)
!==============================================================================!
!Normalization becomes more significant in higher dimension, without a well-
!defined boundry (well) the grid could become increasingly large.
!commented out qMC so we don't need sobol.f90 for convenience
!==============================================================================!
!norm_method          ==>Integration method to normalize p
!xmin(d)              ==>Minimum of normalization domain
!xmax(d)              ==>Maximum of normalization domain
!Nevals               ==>Total number of evaluations
!==============================================================================!
implicit none
character(len=20)::norm_method
integer::Nevals
double precision::xmin(d),xmax(d)
!==============================================================================!
if(norm_method=='uniform_grid'.or.norm_method=='UNIFORM_GRID')then
  call uniform_grid(xmin,xmax,Nevals)
elseif(norm_method=='monte_carlo'.or.norm_method=='MONTE_CARLO')then
  call monte_carlo(xmin,xmax,Nevals)
!elseif(norm_method=='quasi_mc'.or.norm_method=='QUASI_MC')then
!  call quasi_MC(xmin,xmax,Nevals)
else
  stop 'Cannot Identify Normalization Method, Check "normalize_P" Subroutine'
endif
end subroutine normalize_P
!==============================================================================!
subroutine uniform_grid(xmin,xmax,N_1D)
!==============================================================================!
!Compute Integral P with a uniform (square) grid (scales expotentially)
!P(x)~Area_Square/N sum_n=1,N P(x_n)
!==============================================================================!
!xmin(d)              ==>Minimum of normalization domain
!xmax(d)              ==>Maximum of normalization domain
!N_1D                 ==>Number of evaluations along a single dimension
!Ntotal               ==>Total number of evaluations (for all dimensions)
!Moment               ==>First Moment for the distribution
!==============================================================================!
integer::N_1D,Ntotal,i,j
double precision::index1(d),delx(d),xmin(d),xmax(d),x_i(d)
double precision::moment,dummy
!==============================================================================!
!open(20,File='direct_grid.dat')
Moment=0.
Ntotal=(N_1D+1)**d
index1=0
delx(:)=(xmax(:)-xmin(:))/N_1D
do i=1,Ntotal
  do j=1,d
    if(index1(j).eq.N_1D) then
      index1(j)=0
    else
      index1(j)=index1(j)+1
      exit
    endif
  enddo
  x_i(:)=xmin(:)+index1(:)*delx(:)
  Moment=Moment+P_i(x_i)
!  if(V.lt.E_cut) write(20,*) x_i
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
!close(20)
end subroutine uniform_grid
!==============================================================================!
subroutine monte_carlo(xmin,xmax,N_evals)
!==============================================================================!
!Compute Integral P using Monte Carlo (pseudo-random numbers)
!==============================================================================!
!xmin(d)              ==>Minimum of normalization domain
!xmax(d)              ==>Maximum of normalization domain
!Nevals               ==>Total number of evaluations to approximate the integral
!Moment               ==>First Moment for the distribution
!==============================================================================!
implicit none
integer::N_evals,i
double precision::xmin(d),xmax(d),x_trial(d)
double precision::moment,dummy
!==============================================================================!
moment=0d0               !integral_P must be defined to call P, use new variable
do i=1,N_evals
  call random_number(x_trial(:))
  x_trial(:)=xmin(:)+x_trial(:)*(xmax(:)-xmin(:))
  moment=moment+P_i(x_trial)
enddo
dummy=1./N_evals
do i=1,d
  dummy=dummy*(xmax(i)-xmin(i))
enddo
integral_P=moment
integral_P=integral_P*dummy
!==============================================================================!
end subroutine monte_carlo
!==============================================================================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!commenting out sobol routines so I don't need sobol.f90 to compile, these all 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!work still just remove comments if needed
!!!subroutine sobol_unif(skip,x_unif,xmin,xmax)
!!!!==============================================================================!
!!!!Generates a uniform Quasi-Random point spanning rmin-max
!!!!Uses the Sobol Sequence from sobol.f90 made available by John Burkardt
!!!!https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html under GNU LGPL
!!!!==============================================================================!
!!!!skip                 ==>Seed for the random number generator
!!!!x_unif(d)            ==>Uniform Quasi-Random point
!!!!xmin(d)              ==>Minimum of normalization domain
!!!!xmax(d)              ==>Maximum of normalization domain
!!!!==============================================================================!
!!!use sobol
!!!implicit none
!!!integer(kind=8)::skip
!!!double precision::x_unif(d),xmin(d),xmax(d)
!!!!==============================================================================!
!!!x_unif=i8_sobol(int(d, 8), skip)
!!!x_unif=xmin+x_unif*(xmax-xmin)        !Scale uniform distribution to span domain
!!!end subroutine sobol_unif
!!!!==============================================================================!
!!!subroutine quasi_MC(xmin,xmax,N_evals)
!!!!==============================================================================!
!!!!Compute Integral P using quasi-Monte Carlo (quasi-random numbers;Sobol Seq.)
!!!!==============================================================================!
!!!!skip                 ==>Seed for the random number generator
!!!!xmin(d)              ==>Minimum of normalization domain
!!!!xmax(d)              ==>Maximum of normalization domain
!!!!Nevals               ==>Total number of evaluations to approximate the integral
!!!!Moment               ==>First Moment for the distribution
!!!!==============================================================================!
!!!implicit none
!!!integer::N_evals,i
!!!integer(kind=8)::skip
!!!double precision::xmin(d),xmax(d),x_trial(d)
!!!double precision::moment,dummy
!!!!==============================================================================!
!!!moment=0d0               !integral_P must be defined to call P, use new variable
!!!skip=N_evals
!!!do i=1,N_evals
!!!  call sobol_unif(skip,x_trial(:),xmin,xmax)
!!!  moment=moment+P_i(x_trial)
!!!enddo
!!!dummy=1./N_evals
!!!do i=1,d
!!!  dummy=dummy*(xmax(i)-xmin(i))
!!!enddo
!!!integral_P=moment
!!!integral_P=integral_P*dummy
!!!end subroutine quasi_MC
!==============================================================================!
end module morse_mod
!==============================================================================!
program main_grid
use morse_mod
!==============================================================================!
implicit none
character(len=20)::norm_method
integer::N_MMC_box,N_evals,N_MMC_grid,MMC_freq,accept,i,j,k
double precision::time1,time2,Delta_E,dummy,mv_cutoff,Pmax,P_trial
double precision::sigma_trial,E_total
double precision,allocatable,dimension(:)::delx,index1,x_trial,s,x_i
double precision,allocatable,dimension(:)::sigma,x0,xmin,xmax
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
allocate(omega(d))
read(*,*) omega
read(*,*) Npoints
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) N_MMC_box
read(*,*) N_evals
read(*,*) N_MMC_grid
read(*,*) MMC_freq
read(*,*) norm_method
read(*,*) gammas
!==============================================================================!
allocate(x(d,Npoints),sigma(Npoints),x_i(d),delx(d),index1(d))
allocate(x_trial(d),s(d),x0(d),xmin(d),xmax(d))
!==============================================================================!
!                     Box Size for normalizing P (MMC)
!==============================================================================!
call domain_size_P(xmin,xmax,N_MMC_box)
call normalize_P(norm_method,xmin,xmax,N_evals)      !for uniform call with N_1D
!write(*,*) 'Computed normalization'
!==============================================================================!
!                   Generate initial distribution
!==============================================================================!
i=1
x0=0d0
Pmax=P_i(x0)
do while(i.le.Npoints)
  call random_number(x_i)
  x_i(:)=xmin(:)+x_i(:)*(xmax(:)-xmin(:))
  call random_number(dummy)
  if(P_i(x_i)/Pmax.gt.dummy)then                         !Rejection Initial Grid
!  if(V(x_i).lt.E_cut)then                                 !Uniform Initial Grid
     x(:,i)=x_i(:)
     sigma(i)=c_LJ*(P_i(x(:,i))*Npoints)**(-1./d)
     i=i+1
  endif
enddo
E_total=0d0
do i=2,Npoints
   do j=1,i-1
!      write(*,*) Pair_rpl(x(:,i),x(:,j),sigma(i),sigma(j))
      E_total=E_total+Pair_rpl(x(:,i),x(:,j),sigma(i),sigma(j))
   enddo
enddo
open(21,File='grid_ini.dat')
do i=1,Npoints
  write(21,*) x(:,i), P_i(x(:,i)),sigma(i)
enddo
close(21)
!write(*,*) 'Generated random grid'
!==============================================================================!
!                         Generate QRG (greedy search)
!==============================================================================!
open(72,File='grid_conv.dat')
accept=0
mv_cutoff=1.
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  x_trial=x(:,k)+mv_cutoff*sigma(k)*(2*s-1)       !random number: (0,1)==>(-1,1)
  P_trial=P_i(x_trial)
  if(P_trial>1d-20) then                         !Only consider if V(trial)<Ecut
    sigma_trial=c_LJ*(P_trial*Npoints)**(-1./d)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        Delta_E=Delta_E+Pair_rpl(x(:,j),x_trial,sigma(j),sigma_trial) &
        -Pair_rpl(x(:,j),x(:,k),sigma(j),sigma(k))     !Energy change trial move
      endif
    enddo
    if(Delta_E.le.0d0)then
       x(:,k)=x_trial(:)
       sigma(k)=sigma_trial
       accept=accept+1
       E_total=E_total+Delta_E
    endif
  endif
!for MMC want acceptance ~30-50%, adjust trial movement displacement accordingly
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/MMC_freq.lt.0.4)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
    accept=0
    write(72,*) i,mv_cutoff,E_total
  endif
enddo
close(72)
open(22,File='grid.dat')
do i=1,Npoints
  write(22,*) x(:,i)
enddo
close(22)
!==============================================================================!
open(98,file='intp.dat')
write(98,*) N_evals, integral_P
close(98)
!==============================================================================!
!                                 Output
!==============================================================================!
call cpu_time(time2)
open(99,file='out.txt')
write(99,*) 'd ==> ', d
write(99,*) 'Pmax=', Pmax
write(99,*) 'E_total=',E_total
do i=1,d
  write(99,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_evals ==> ',N_evals
write(99,*) 'Npoints ==> ',Npoints
write(99,*) 'c_LJ ==> ',c_LJ
write(99,*) 'N_MMC_box ==> ',N_MMC_box
write(99,*) 'N_MMC_grid ==> ',N_MMC_grid
write(99,*) 'MMC_freq ==> ',MMC_freq
write(99,*) 'gammas==> ',gammas
write(99,*) 'Simulation Time==> ',time2-time1
close(99)
end program main_grid
