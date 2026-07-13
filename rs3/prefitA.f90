! prefitA.f90
! module for initial estimation of max_blaze polynomial for resimplex3 program
! the estimation will be based on results of the presimplex module
!
! Pavel Dolezal October 6th 2025
! pavel.dolezal@matfyz.cuni.cz

module prefitA_module

contains

subroutine prefitA()

use variables_module
use loader_module
use prefitA_func_module
use amoeba_mod

implicit none

integer :: i,j,iter, n, posmin, ndim, k, l
double precision, allocatable :: x(:), p(:,:), e(:), z(:), chi2

print *, "prefitA started"
flush(6)

ndim = max_ord+1

allocate (x(ndim))
allocate (e(ndim))
allocate (p(ndim+1, ndim))
allocate (z(ndim+1))

do i=1, n_dat
	e(:) = emb
	dataset = i

	! initialization of x
	x(:) = mb(part,i,:)
	
	do l=1,n_conv
		print *
		print *, "prefitA, part, dataset, convergence: ", part, i ,l
		print *, "x = ", x
		print *, "e = ", e
		flush(6)

		! create initial simplex
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
			z(k) = prefitA_func(x)
			print *, "z(", k, ") = ", z(k)
		end do

		r = 0 ! reset counter for chi2_func
		n = 0 ! reset counter of non skipped orders for chi2_func	
		iter = 0

		! let's roll
		call amoeba(p,z,ftol_pres,prefitA_func,iter,200000)

		! determine the best result
		posmin = minloc(z,1)
		x = p(posmin,:)
		chi2 = z(posmin)
		print *, "prefitA, current dataset, part, convergence, chi2: ", i, part, l,chi2
		print *, "x = ", x

		! save results for this dataset
		mb(part,i,:) = x(:)
		e=e*0.8
	end do ! l=n_conv

	! fake rectification - call prefitA_func to check resulted fit (fit_checkA.plt) 
        rectifinal = .true.
        chi2 = prefitA_func(x)
        rectifinal = .false.

end do ! i=n_dat

end subroutine prefitA

end module prefitA_module
