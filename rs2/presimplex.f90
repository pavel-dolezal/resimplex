! presimplex.f90
! module for finding lambda_c and max_blaze for each order and each dataset
! lambda_c and max_blaze are free, meaning they are constant for each order and dataset
! and are not bind by any function
! Pavel Dolezal February 19th 2025
! pavel.dolezal@matfyz.cuni.cz

module presimplex_module

contains

subroutine presimplex

use variables_module
use loader_module
use chi2_func_presx_module
use amoeba_mod

implicit none

integer :: iter, n, posmin, ndimr, i,j,k
double precision :: xr(2), pr(3,2), er(2), zr(3), chi2

print *, "resimplex started"
flush(6)

ndimr = 2

DO i=1, n_dat
	cur_dat = i
	do j=1, n_ord	
	        if ((skipped_orders(j) == 0).and.((j>min2).and.(j<=max2))) then
			cur_ord=j		
			order = j !integer -> real, for chi2_func_presx	

			! initialization of x
			xr(1) = lc(i,j)
			xr(2) = m(i,j)
			print *, xr
			er(1) = elc
			if (j>=n_ord-15) then
				er(2) = em/10
			else
				er(2) = em
			end if

			print *, "er = ",  er
			flush(6)

			! create initial simplex
			do k=1, ndimr+1
				pr(k,:)=xr
				if (k<ndimr+1) pr(k,k)=pr(k,k)+er(k)
			end do
		
			print *, "The initial simplex:"
			do k=1, ndimr+1  
				write(*, *) "pr(",k,",:) = ", pr(k,:)
			end do
			flush(6)
			! determine chi2 of the initial simplex
			do k=1, ndimr+1
				xr=pr(k,:)
				zr(k) = chi2_func_presx(xr)
				print *, "zr(", k, ") = ", zr(k)
			end do

			r = 0 ! reset counter for chi2_func_presx
			n = 0 ! reset counter of non skipped orders for chi2_func_presx	
			iter = 0

			! let's roll
			call amoeba(pr,zr,ftol_pres,chi2_func_presx,iter,ITMAX1)

			! determine the best result
			posmin = minloc(zr,1)
			xr = pr(posmin,:)
			chi2 = zr(posmin)
			print *, "presimplex, current dataset, current dataset, chi2: ",i, j, chi2
			print *, "xr = ", xr
			lc(i,j) = xr(1)
			m(i,j) = xr(2)
		end if ! skipping orders
	end do ! j=n_ord
END DO ! i=n_dat

end subroutine presimplex

end module presimplex_module
