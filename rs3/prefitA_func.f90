! prefitA_func.f90
! module for bsimplex program

! Pavel Dolezal August 14th 2025
! pavel.dolezal@matfyz.cuni.cz
module	prefitA_func_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function prefitA_func(x)

use loader_module
!use blaze_module
use variables_module

implicit none
double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, exponent

double precision, parameter :: pi = 4.0d0*Datan(1.0d0)
integer :: i, j, k, o, ii

! chi2 computation
chi2 = 0.d0
do j=1, n_ord
        if ((skipped_orders(j) == 0).and.((j>min2).and.(j<=max2))) then
	!if ((j>min2).and.(j<=max2)) then
		i=dataset
		order=j

		max_blazef = 0.d0
		do ii=1, max_ord+1
			exponent = ii-1
			max_blazef = max_blazef + x(ii)*order**(exponent)
		end do
		
		! write down the resulted fit to fit_checkA.asc
		if (rectifinal) then
			write (667,*) j, free_mb(i,j), max_blazef, i
		end if

		chi2 = chi2 + (free_mb(i,j)-max_blazef)**2
	end if
end do ! j=n_ord

! printing every 10th iteration to chi2_func.tmp
if (tmp_prefitA_func) then
	if (r==10) then
		call cpu_time(t2)
		write(56, *) x, chi2, t2-t1
		r=0
	else
		r=r+1
	end if
end if

prefitA_func = chi2

return

end function prefitA_func

end module prefitA_func_module
