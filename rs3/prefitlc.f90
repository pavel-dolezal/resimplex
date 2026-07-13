! prefitlc.f90
! Pavel Dolezal August 14th 2025
! pavel.dolezal@matfyz.cuni.cz

module prefitlc_module

contains

subroutine prefitlc()

use variables_module
use loader_module
use prefitlc_func_module
use amoeba_mod

implicit none

integer :: i,j,iter, n, posmin, ndim, k, l
double precision, allocatable :: x(:), p(:,:), e(:), z(:), chi2

print *, "prefitlc started"
flush(6)

ndim = lc_ord+1

allocate (x(ndim))
allocate (e(ndim))
allocate (p(ndim+1, ndim))
allocate (z(ndim+1))

do i=1, n_dat
	e = elc
	dataset = i

	! initialization of x
	x(:) = lc(part,i,:)
	
	do l=1,n_conv
		! create initial simplex

		print *	
		print *, "prefitlc, part, dataset, convergence: ", part, i, l
		print *, "x = ", x
		print *, "e = ", e
		flush(6)

		do k=1, ndim+1
			p(k,:)=x
			if (k<ndim+1) p(k,k)=p(k,k)+e(k)
		end do

		print *, "The initial simplex:"
		do k=1, ndim+1  
			write(*, *) "p(",k,",:) = ", p(k,:)
		end do
		flush(6)

		! determine chi2 of the initial simplex
		do k=1, ndim+1
			x=p(k,:)
			z(k) = prefitlc_func(x)
			print *, "z(", k, ") = ", z(k)
		end do

		r = 0 ! reset counter for chi2_func
		n = 0 ! reset counter of non skipped orders for chi2_func	
		iter = 0

		! let's roll
		call amoeba(p,z,ftol_pres,prefitlc_func,iter,200000)

		! determine the best result
		posmin = minloc(z,1)
		x = p(posmin,:)
		chi2 = z(posmin)
		print *, "prefitlc, current dataset, part, convergence, chi2: ", i, part, l,chi2
		print *, "x = ", x

		! save results for this dataset
		lc(part,i,:) = x(:)	
	end do

	! false rectification - check resulted fit (fit_check.plt)
        rectifinal = .true.
        chi2 = prefitlc_func(x)
        rectifinal = .false.
	
end do ! i=n_dat

close(50)

end subroutine prefitlc

end module prefitlc_module
