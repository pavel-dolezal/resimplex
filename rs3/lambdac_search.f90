! lambdac_search.f90
! module for finding lambda_c and max_blaze for each order and each dataset
! lambda_c and max_blaze are free, meaning they are constant for each order and dataset
! and are not bind by any function
!
! Pavel Dolezal February 19th 2025
! pavel.dolezal@matfyz.cuni.cz

module lambdac_search_module

contains

subroutine lambdac_search(min,max)

use variables_module
use loader_module
use chi2_funcA_module
use amoeba_mod
use prefitlc_module
use prefitA_module

implicit none

!double precision, intent(inout) :: lc(:,:)
integer :: i,j,iter, n, posmin, ndim, k, it
double precision, allocatable :: x(:), p(:,:), e(:), z(:), chi2, chi2A, chi2B
integer, intent(IN) :: min, max

print *, "lambdac_search started"
flush(6)

min2=min
max2=max

ndim = size(lc,3)

allocate (x(ndim))
allocate (e(ndim))
allocate (p(ndim+1, ndim))
allocate (z(ndim+1))

! fit free max_blaze and lambda_c with polynomials
call prefitlc
call prefitA

chi2_all_ls = 0.d0
do i=1, n_dat
	e(:) = elc
	dataset = i
	
        if (rectify) then
                write(out_name, "(A,A,I5.5,A)") "lambdac_check/", trim(prefix), list(i), ".asc"
                print *, out_name, "opened"
		if (part==1) then
	                open(unit = 300, file = out_name, status = "replace")
			write (300,*) "# wavelength [A] & normalised_intensity I_lambda [] & sigma [] ",&
	                "& interpolated? [logical] & blaze_value [] & order_parity []"
	               write (300,*) "#",comment
		else
			open(unit = 300, file = out_name, status = "old",access = 'append')
		end if

        end if

	! initialization of x
	x(:) = lc(part,i,:)
	chi2A=0.d0
	chi2B=0.d0
	it = 0
	do !until chi2A-chi2B < conv_treshold
		print *
		it = it+1
		print *, "start of lambdac_search: part, dataset,convergence = ", part, i, it
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
			z(k) = chi2_funcA(x)
			print *, "z(", k, ") = ", z(k)
		end do

		r = 0 ! reset counter for chi2_func
		n = 0 ! reset counter of non skipped orders for chi2_func	
		iter = 0
		
		print *, "lambdac_search: part, dataset", part, i

		! let's roll
		call amoeba(p,z,ftol_pres,chi2_funcA,iter,ITMAX2)

		! determine the best result
		posmin = minloc(z,1)
		x = p(posmin,:)
		chi2 = z(posmin)
		
		print *, "x = ", x
		
		! stop the loop if final convergence has reduced chi2 by less than conv_treshold
                chi2B = chi2
                
		if ((chi2A>0).and.((chi2A-chi2B)<conv_treshold)) then
                        print *, "the final convergence #", it
                        print *, "chi2A-chi2B = ", chi2A-chi2B
                        exit
                else
                        chi2_diff = chi2A-chi2B
                        chi2A = chi2B
                end if

		if (it>5) then
			e=e*0.8
		end if		

                print *, "part ",part, ", convergence: ", it, ", last chi2A-chi2B = ", chi2_diff, "dataset = ", dataset

	end do
	! save results for this dataset
	lc(part,i,:) = x(:)
	
	! call chi2_func for final rectification
	if (rectify) then
	        rectifinal = .true.
	        chi2 = chi2_funcA(x)
	        rectifinal = .false.
		chi2_all_ls = chi2_all_ls + chi2_allA
	end if

	print *, "end of lambdac_search: part, current dataset, chi2: ", part, i, chi2

end do ! i=n_dat

close(300)

deallocate (x)
deallocate (e)
deallocate (p)
deallocate (z)

end subroutine lambdac_search

end module lambdac_search_module
