!=============================================================================80
!                         Histogram of Gaussian Grid
!=============================================================================80
!       Discussion:
!Generate a histogram for a distribution function: P(r) (gaussain example)
!Requires a set of points to
!==============================================================================!
!       Modified:                                                              !
!   14 July 2021                                                               !
!       Author:                                                                !
!   Shane Flynn                                                                !
!==============================================================================!
module hist_mod
implicit none
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
!NG              ==>Number of points to bin
!d               ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!rrange          ==>(d) domain range
!==============================================================================!
integer,  parameter :: d = 2         !should work for any d, tested for d=2 only
integer :: NG
double precision :: rrange(d)
!==============================================================================!
contains
!==============================================================================!
function P(r)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(r)
!r              ==>(d) ith points coordinate r^i_1,..,r^i_d
!==============================================================================!
implicit none
double precision :: r(d),P
double precision, parameter :: pie = acos(-1d0)
!==============================================================================!
P = (2.*pie)**(-d/2.)*exp(-0.5*sum(r(:)**2))
end function P
!==============================================================================!
subroutine find_range(r)
!==============================================================================!
!Domain Range
!==============================================================================!
!rmin           ==>(d) minimum of the range domain
!rmax           ==>(d) maximum of the range domain
!r              ==>(d,NG) ith points coordinate r^i_1,..,r^i_d
!==============================================================================!
implicit none
integer :: i, k
double precision :: rmin(d), rmax(d), r(d,NG)
!==============================================================================!
rmin(:) = r(:,1)
rmax(:) = r(:,1)
do i=2,NG
   do k=1,d
      if(rmin(k) > r(k,i)) rmin(k) = r(k,i)
      if(rmax(k) < r(k,i)) rmax(k) = r(k,i)
   enddo
enddo
do k=1,d
   write(*,*) k,'  range=',rrange(k)
   rrange(k) = 1/(rmax(k)-rmin(k))
enddo
end subroutine find_range
!==============================================================================!
subroutine histogram(r)
!==============================================================================!
!Generate a histogram for a given set of points
!==============================================================================!
implicit none
integer :: i, j, k
double precision, parameter :: rmin=0d0, rmar=2d0
double precision, parameter :: dr=(rmar-rmin)/50
double precision :: hist(1:100), r(d,NG), dist_aver, dist, dist1, Pr
do i=1,NG
   Pr = P(r(:,i))
   dist1 = 1d10
   do j=1,NG
      if(i.ne.j) then
         dist = sqrt(sum((rrange(:)*(r(:,i)-r(:,j)))**2))
         if(dist < dist1) dist1=dist
      endif
   enddo
   dist_aver = dist_aver+dist1*Pr**(1./d)
enddo
dist_aver = dist_aver/NG
write(*,*) 'dist_aver=', dist_aver
hist=0d0
do i=1,NG
   Pr = P(r(:,i))
   do j=1,NG
      if(i.ne.j) then
         dist=sqrt(sum((rrange(:)*(r(:,i)-r(:,j)))**2))*Pr**(1./d)/dist_aver
         k = (dist-rmin)/dr+1
         if(k>100) k=100
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
end module hist_mod
!==============================================================================!
program hist
use  hist_mod
!==============================================================================!
implicit none
integer :: i
double precision,allocatable :: r(:,:)
!==============================================================================!
!                              Read Input File                                 !
!==============================================================================!
read(*,*) NG
 !==============================================================================!
allocate(r(d,NG))
!==============================================================================!
!                               read in grid
!==============================================================================!
open(21,file='grid.dat')
do i=1,NG
   read(21,*) r(:,i)
enddo
close(21)
!==============================================================================!
!                         Generate QRG (greedy search)
!==============================================================================!
call find_range(r)
write(*,*) 'Compute the pair correlation function'
open (8, file='hist.dat')
call histogram(r)
close(8)
end program hist
