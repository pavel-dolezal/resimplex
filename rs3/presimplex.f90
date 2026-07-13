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
use chi2_func2_module
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
	
		if (j<split(1)) then
			part = 1
		else
			if (j<split(2)) then
				part = 2	
			else
				if (j<split(3)) then
					part = 3
				else
					part = 4
			
				end if
			end if
		end if
		
		cur_ord=j		
		order = j !integer -> real, for chi2_func2	
		print *	
		print *, "presimplex, current dataset, current order, chi2: ",i, j, chi2

		! initialization of x
		xr(1) = free_lc(i,j)
		xr(2) = free_mb(i,j)
		print *, "initial xr = ", xr
		er(1) = elc_presx
		if (j>=n_ord-15) then
			er(2) = em_presx/10
		else
			er(2) = em_presx
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
			zr(k) = chi2_func2(xr)
			print *, "zr(", k, ") = ", zr(k)
		end do

		r = 0 ! reset counter for chi2_func2
		n = 0 ! reset counter of non skipped orders for chi2_func2	
		iter = 0

		! let's roll
		call amoeba(pr,zr,ftol_pres,chi2_func2,iter,ITMAX1)

		! determine the best result
		posmin = minloc(zr,1)
		xr = pr(posmin,:)
		chi2 = zr(posmin)
		print *, "final xr = ", xr
		print *, "chi2 = ", chi2
		free_lc(i,j) = xr(1)
		free_mb(i,j) = xr(2)
	end do ! j=n_ord
END DO ! i=n_dat

end subroutine presimplex

end module presimplex_module
