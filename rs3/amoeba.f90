! amoeba in F90 style
module amoeba_mod

contains

subroutine amoeba(p,y,ftol,func,iter,ITMAX)

use nrtype; use nrutil, only : assert_eq,imaxloc,iminloc,nrerror,swap

implicit none
integer(i4b), intent(out) :: iter
integer(I4B), intent(in) :: ITMAX
real(DP), intent(in) :: ftol
real(DP), dimension(:), intent(inout) :: y
real(DP), dimension(:,:), intent(inout) :: p

interface

function func(x)
use nrtype
use loader_module
use variables_module
implicit none
real(DP), dimension(:), intent(in) :: x
real(DP) :: func
end function func
end interface

real(DP), parameter :: TINY=1.0d-10
! Minimization of the function func in N dimensions by the downhill simplex method of
! Nelder and Mead. The (N + 1) × N matrix p is input. Its N + 1 rows are N -dimensional
! vectors that are the vertices of the starting simplex. Also input is the vector y of length
! N + 1, whose components must be preinitialized to the values of func evaluated at the
! N + 1 vertices (rows) of p; and ftol the fractional convergence tolerance to be achieved
! in the function value (n.b.!). On output, p and y will have been reset to N + 1 new points
! all within ftol of a minimum function value, and iter gives the number of function
! evaluations taken.
! Parameters: The maximum allowed number of function evaluations, and a small number.
integer(I4B) :: ihi,ndim
! Global variables.
real(DP), dimension(size(p,2)) :: psum
call amoeba_private

contains

subroutine amoeba_private

implicit none
integer(I4B) :: i,ilo,inhi
real(DP) :: rtol,ysave,ytry,ytmp

ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
iter=0
psum(:)=sum(p(:,:),dim=1)
do ! Iteration loop.
  ilo=iminloc(y(:)) ! Determine which point is the highest (worst),
  ihi=imaxloc(y(:)) ! next-highest, and lowest (best).
  ytmp=y(ihi)
  y(ihi)=y(ilo)
  inhi=imaxloc(y(:))
  y(ihi)=ytmp
  rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
  ! Compute the fractional range from highest to lowest and return if satisfactory.
  if (rtol < ftol) then ! If returning, put best point and value in slot 1.
    call swap(y(1),y(ilo))
    call swap(p(1,:),p(ilo,:))
    return
  end if
  if (iter >= ITMAX) then
!    call nrerror('ITMAX exceeded in amoeba')
    write(*,*) '# Warning: ITMAX exceeded in amoeba'
    return
  endif
! Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex
! across from the high point, i.e., reflect the simplex from the high point.
  ytry=amotry(-1.0_dp)
  iter=iter+1
  if (ytry <= y(ilo)) then ! Gives a result better than the best point, so
    ytry=amotry(2.0_dp)    ! try an additional extrapolation by a factor of 2.
    iter=iter+1
  else if (ytry >= y(inhi)) then ! The reflected point is worse than the sec-
    ysave=y(ihi)                 ! ond highest, so look for an intermediate
    ytry=amotry(0.5_dp)          ! lower point, i.e., do a one-dimensional contraction.
    iter=iter+1
    if (ytry >= ysave) then  ! Can’t seem to get rid of that high point. Better contract around the lowest (best) point.
      p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
      do i=1,ndim+1
        if (i /= ilo) y(i)=func(p(i,:))
      end do
      iter=iter+ndim  ! Keep track of function evaluations.
      psum(:)=sum(p(:,:),dim=1)
    end if
  end if
end do  ! Go back for the test of doneness and the next iteration.

end subroutine amoeba_private

function amotry(fac)

implicit none
real(DP), intent(in) :: fac
real(DP) :: amotry
! Extrapolates by a factor fac through the face of the simplex across from the high point,
! tries it, and replaces the high point if the new point is better.
real(DP) :: fac1,fac2,ytry
real(DP), dimension(size(p,2)) :: ptry
fac1=(1.0_dp-fac)/ndim
fac2=fac1-fac
ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
ytry=func(ptry) ! Evaluate the function at the trial point.
if (ytry < y(ihi)) then ! If it’s better than the highest, then replace
  y(ihi)=ytry           ! the highest.
  psum(:)=psum(:)-p(ihi,:)+ptry(:)
  p(ihi,:)=ptry(:)
end if
amotry=ytry

end function amotry

end subroutine amoeba

end module amoeba_mod
