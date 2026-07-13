! prefit_t.f90
! module for resimplex program
! this module fit one of resulting values of adetg blaze parameters from resimplex1
! this fit is susequently used as first estimation of adetg polynomials in the main program
!
! Pavel Dolezal August 14th 2025
! pavel.dolezal@matfyz.cuni.cz

module prefit_t_module

contains

subroutine prefit_t()

use variables_module
use loader_module
use prefit_t_func_module
use amoeba_mod

implicit none

integer :: i,j,iter, n, posmin, ndim2, k, l, n_conv, part2
double precision, allocatable :: x(:), p(:,:), e(:), z(:), chi2

print *, "prefit_t started"
flush(6)

ndim2 = (adtge_poly+1)

allocate (x(ndim2))
allocate (e(ndim2))
allocate (p(ndim2+1, ndim2))
allocate (z(ndim2+1))

e(:) = et

print *, "e = ", e
n_conv = 6

do part2=1, 1
	
	if (part2==1) then	
		fit_min = 0
	else
		fit_min = split(part2-1)
	end if 

	if (part2 == 4) then
		fit_max = n_ord+1
	else
		fit_max = split(part2)
	end if

	! initialization of x
	x(:) = t(part2,:)

	print *, "prefit_t, part2 =  ", part2
	print *, "x = ", x
	print *, "e = ", e
	flush(6)
	do l=1,n_conv
		! create initial simplex

		do k=1, ndim2+1
			p(k,:)=x
			if (k<ndim2+1) p(k,k)=p(k,k)+e(k)
		end do

		print *, "The initial simplex:"
		do k=1, ndim2+1  
			write(*, *) "p(",k,",:) = ", p(k,:)
		end do
		flush(6)

		! determine chi2 of the initial simplex
		do k=1, ndim2+1
			x=p(k,:)
			z(k) = prefit_t_func(x)
			print *, "z(", k, ") = ", z(k)
		end do

		r = 0 ! reset counter for chi2_func
		n = 0 ! reset counter of non skipped orders for chi2_func	
		iter = 0

		! let's roll
		call amoeba(p,z,ftol_pres,prefit_t_func,iter,200000)

		! determine the best result
		posmin = minloc(z,1)
		x = p(posmin,:)
		chi2 = z(posmin)
		print *, "prefit_t, current part2, convergence, chi2: ", part2,l,chi2
		print *, "x = ", x

		! save results for this dataset
		do i=1, 4
			t(i,:) = x(:)
		end do
	end do

	write(135, *) "t = ", t(1,:)

	! call chi2_func for final rectification
	if (rectify) then
	        rectifinal = .true.
	        chi2 = prefit_t_func(x)
	        rectifinal = .false.
	end if
end do

end subroutine prefit_t

end module prefit_t_module
