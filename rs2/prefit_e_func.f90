! prefit_e_func.f90
! module for prefit module of the resimplex program
!
! Pavel Dolezal October 13th 2025
! pavel.dolezal@matfyz.cuni.cz

module	prefit_e_func_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function prefit_e_func(x)

use loader_module
use variables_module

implicit none
double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, exponent

double precision, parameter :: pi = 4.0d0*Datan(1.0d0)
integer :: i, j, k, o, ii, jj

! chi2 computation
chi2 = 0.d0

do j=1, n_ord
        if ((skipped_orders(j) == 0).or.(rectifinal)) then
		order=j

		epsilonf = 0.d0
                do i=1, adtge_poly+1
                        exponent = i-1
                        epsilonf = epsilonf+x(i)*order**(exponent)
                end do

		chi2 = chi2 + (free_eps(j)-epsilonf)**2

		if (rectifinal) then
			write (668,*) j, free_eps(j),epsilonf
		end if
	end if
end do ! j=n_ord

! printing current iteration to chi2_func_fit.tmp
if (tmp_prefit_func) then
	if(r==10) then
		call cpu_time(t2)
		write(57, *) x, chi2, t2-t1
		r=0
	else
		r=r+1
	end if
end if

prefit_e_func = chi2
return

end function prefit_e_func

end module prefit_e_func_module
