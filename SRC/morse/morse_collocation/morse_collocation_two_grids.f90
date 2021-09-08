!=============================================================================80
!                  ND-Morse EigenSpectra (Sparse Quadrature)
!==============================================================================!
!Compute the EigenSpectra for an nD-Morse Potential
!Quasi-Regular Grid as Input; generate (Gaussian Widths) from nearest neghbor
!This code uses LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   30 November 2020
!==============================================================================!

module QRGC_mod
  implicit none
  !==============================================================================!
  !                            Global Variables
  !==============================================================================!
  !d              ==>i-th gaussian dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
  !NG             ==>Number of Gaussian to generate
  !omega          ==>(d) Parameter for Morse Potential
  !==============================================================================!
  integer::d,NG,Nc
  double precision, parameter :: pi=4.*atan(1d0)
  double precision, allocatable, dimension(:) :: omega,alpha
  double precision, allocatable :: x(:,:), r(:,:)
  !==============================================================================!
contains

  function Phi(i,x)
    ! Normalized d-dimensional Gaussian
    implicit none
    integer :: i
    double precision :: x(d),Phi
    Phi=exp(-alpha(i)*sum((r(:,i)-x(:))**2))   !removed normalization
!    Phi=(2*alpha(i)/pi)**d/4.0*exp(-alpha(i)*sum((r(:,i)-x(:))**2))   !remove normalization later
  end function Phi
    
  !==============================================================================!
  function V(x)
    !==============================================================================!
    !Sum of 1D-Morse Potential
    !==============================================================================!
    !x              ==>(d) ith atoms coordinates
    !V              ==>evaluate V(x_i)
    !D_morse        ==>Parameter for Morse Potential
    !==============================================================================!
    implicit none
    double precision::x(d),V
    double precision,parameter::D_morse=12.
    V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
  end function V

  function T(i,n)
    !Using 5-point finite difference to compute the laplassian of Phi
    double precision,parameter::ss=1d-5                               !need to check
    double precision, parameter :: a(-2:2)=(/1.,-16.,30.,-16.,1./)/(24*ss**2)
    integer :: i,k,l,n,m
    double precision :: T, y(d)
    T=0d0
    y(:)=x(:,n)
    do k=1,d
       do l=-2,2     
          y(k)=x(k,n)+l*ss
          T=T+a(l)*Phi(i,y)
          y(k)=x(k,n)
       enddo
    enddo
  end function T



  subroutine SF_to_int(x,q,flag)
    !   convert Space-fixed (x) into internal (q) coordinates
    !   In this special case the dimensions are the same and the transformation is the identity
    double precision :: x(d),q(d)
    logical:: flag
    if(flag) then
       q=x
    else
       x=q
    endif
  end subroutine SF_to_int


  !==============================================================================!  
  subroutine matrix_elements_collocation(Smat,Hmat)
    !==============================================================================!
    ! Using two grids
    integer :: i,j,n,k,l
    double precision :: Smat(NG,NG),Hmat(NG,NG),V_(Nc),H_i(Nc),Phi_i(Nc),Phi_jn,y(d)
    !Using 5-point finite difference to compute the laplassian of Phi
    double precision,parameter::ss=1d-5                               !need to check
    double precision, parameter :: a(-2:2)=(/1.,-16.,30.,-16.,1./)/(24*ss**2)
    do n=1,Nc
       V_(n)=V(x(:,n))
    enddo
    Smat=0d0
    Hmat=0d0
    do i=1,NG
       do n=1,Nc
          y(:)=x(:,n)
          Phi_i(n)=Phi(i,y)
          H_i(n)=(V_(n)+d*a(0))*Phi_i(n)
          do k=1,d
             do l=-2,2
                if(l.ne.0) then
                   y(k)=x(k,n)+l*ss
                   H_i(n)=H_i(n)+a(l)*Phi(i,y)
                   y(k)=x(k,n)
                endif
             enddo
          enddo
       enddo
       do j=1,NG
          do n=1,Nc
             Phi_jn=Phi(j,x(:,n))
             Smat(i,j)=Smat(i,j)+Phi_i(n)*Phi_jn
             Hmat(i,j)=Hmat(i,j)+H_i(n)*Phi_jn
          enddo
       enddo
    enddo
    ! Symmetrize
    do i=2,NG
       do j=1,i-1
          Hmat(i,j)=(Hmat(i,j)+Hmat(j,i))/2    ! i>j
          Hmat(j,i)=Hmat(i,j)
       enddo
    enddo
  end subroutine matrix_elements_collocation
  
    subroutine hpsort(N,A,B)
    Implicit None
    Integer, Intent(In) :: N
    Double Precision, Intent(Inout) :: A(N),B(N)
    Integer :: I,IR,J,L
    Double Precision :: RA,RB
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
20  if (J.le.IR) then
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
end module QRGC_mod
      !==============================================================================!
      program QRGC_Spec
        use QRGC_mod
        !==============================================================================!
        !unif_grid      ==>If True set alpha to be a constant, else use nearest neighbor
        !grid_in        ==>Filename Containing Gridpoints
        !theory_in      ==>Filename Containing Analytic Eigenvalues
        !NG             ==>Number of Gaussian Basis Functions (gridpoints)
        !GH_order       ==>Number of Points for evaluating the potential (Gauss-Hermite)
        !alpha0         ==>Flat Scaling Parameter for Gaussian Widths
        !alpha          ==>(d) Gaussian Widths
        !RCN            ==>Recriprical Convergence Number, stability of Overlap Matrix
        !x              ==>(d) ith atoms coordinates
        !x_ij           ==>i-jth gaussian center (product of Gaussians is a Gaussian)
        !Smat           ==>(NG,NG) Overlap Matrix
        !Hmat           ==>(NG,NG) Hamiltonian Matrix
        !==============================================================================!
        implicit none
        character(len=100)::Gaussians_centers,Collocation_points,theory_in
        integer::GH_order,sg_level,quad_size,i,j,k,tsize
        external gqn,gqn_order                                        !sparse_quadrature
        integer,allocatable,dimension(:)::l
        double precision::x2,alpha0,time1,time2,time3
        double precision,allocatable,dimension(:) :: Er,Ei,theory,BETA
        double precision,allocatable,dimension(:,:) :: Smat,Hmat,VL,VR
        !==============================================================================!
        !                            LLAPACK dsygv variables
        !==============================================================================!
        integer::itype,info,lwork
        double precision,allocatable,dimension(:)::work
        !==============================================================================!
        !                             Read Input Data File
        !==============================================================================!
        call cpu_time(time1)
        read(*,*) d
        allocate(omega(d))
        read(*,*) omega
        read(*,*) Gaussians_centers,NG
        read(*,*) Collocation_points,Nc
        Nc=Nc+NG
        read(*,*) theory_in
        read(*,*) tsize     !how many theory eigenvalues to read in
        read(*,*) alpha0

        allocate(x(d,Nc),r(d,NG),alpha(NG),Smat(NG,NG),Hmat(NG,NG))

        !==============================================================================!
        !                           Read GridPoints x(d,NG)
        !==============================================================================!
        open(17,File=Gaussians_centers)
        do i=1,NG
           read(17,*) r(:,i)
           x(:,i)=r(:,i)
        enddo
        close(17)
        open(17,File=Collocation_points)
        do i=NG+1,Nc
           read(17,*) x(:,i)
        enddo
        close(17)
        !==============================================================================!
        !         Symmetric gaussians. Use nearest neighbor to determine alpha
        !==============================================================================!
        do i=1,NG
           alpha(i)=1d20                            !large distance for placeholder
           do j=1,NG
              if(j.ne.i) then
                 x2=sum((r(:,i)-r(:,j))**2)          !distance between gridpoints
                 if(x2<alpha(i)) alpha(i)=x2
              endif
           enddo
           alpha(i)=alpha0/alpha(i)
        enddo
        open(unit=18,file='alphas.dat')
        do i=1,NG
           write(18,*) alpha(i)
        enddo
        close(18)
        call matrix_elements_collocation(Smat,Hmat)
        call cpu_time(time2)
        write(*,*) 'Constructed the matrices time ==>', time2-time1
        deallocate(x,r,alpha)
        !==============================================================================!
        !               Eigenvalues of the Hamiltonian matrix
        !==============================================================================!
        write(*,*)  "Solving the generalized eigenvalue problem"
        open(unit=20,file='eigenvalues.dat')


        info=0
        itype=1
        Lwork=max(1,3*NG-1)
        allocate(work(Lwork),Er(NG))
        call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,Er,work,Lwork,info)
        write(*,*) 'info ==> ', info
        do i=1, NG
           write(20,*) Er(i)
        enddo


!!$        lwork=max(1,8*NG)
!!$        allocate(work(lwork),Er(NG),Ei(NG),BETA(NG))
!!$        itype=1
!!$        call dggev('n','n',NG,Hmat,NG,Smat,NG,Er,Ei,BETA,VL,NG,VR,NG,work,lwork,info)
!!$        deallocate(work,Hmat,Smat)
!!$        Er(:)=Er(:)/BETA(:)
!!$        Ei(:)=Ei(:)/BETA(:)
!!$        call hpsort(NG,Er,Ei)
!!$        write(*,*) 'info ==> ', info
!!$        open(unit=20,file='eigenvalues.dat')
!!$        do i=1, NG
!!$           write(20,*) Er(i)
!!$        enddo

        close(20)
        !==============================================================================!
        !                              Exact Eigenvalues
        !==============================================================================!
        allocate(theory(tsize))
        open(21,File=theory_in)
        do i=1,tsize
           read(21,*) theory(i)
        enddo
        close(21)
        open(unit=22,file='abs_error.dat')
        open(unit=23,file='rel_error.dat')
        open(unit=24,file='alpha_abs_error.dat')
        open(unit=25,file='alpha_rel_error.dat')
        do i=1,tsize
           write(22,*) i, abs(theory(i)-Er(i)),Ei(i)
           write(23,*) i, abs(Er(i)-theory(i))/theory(i),Ei(i)/theory(i)
           write(24,*) alpha0, abs(theory(i)-Er(i))
           write(25,*) alpha0, (Er(i)-theory(i))/theory(i)
        enddo
        close(22)
        close(23)
        close(24)
        close(25)
        open(unit=26,file='alpha_rel_150.dat')
        do i=1,150
           write(26,*) alpha0, (Er(i)-theory(i))/theory(i)
        enddo
        close(26)
        !==============================================================================!
        !                                Output file                                   !
        !==============================================================================!
        call cpu_time(time3)
        write(*,*) 'Diagonalization time ==>', time3-time2
        write(*,*) 'Total time ==>', time3-time1
        open(99,file='simulation.dat')
        write(99,*) 'dimensionality ==> ', d
        write(99,*) 'omega ==> ', omega
        write(99,*) 'NG ==> ', NG
        write(99,*) 'Nc ==> ', NG
        write(99,*) 'alpha0==>', alpha0
        write(99,*) 'Matrix conxtruction time ==>', time2-time1
        write(99,*) 'Diagonalization time ==>', time3-time2
        write(99,*) 'Total time ==>', time3-time1
        close(99)
        write(*,*) 'Hello Universe!'
      end program QRGC_Spec
