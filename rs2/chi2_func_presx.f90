! chi2_func_presx.f90
! module for resimplex program
! computation of chi2 for subroutine presimplex and for given parameters x(...) of the blaze function
! x contains just one lambda_c and one max_blaze for given order and dataset

! Pavel Dolezal September 9th 2025
! pavel.dolezal@matfyz.cuni.cz
module	chi2_func_presx_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function chi2_func_presx(x)

use loader_module
use blaze_module
use variables_module

implicit none

double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, exponent

integer :: i, j, k, ii

double precision, parameter :: pi = 4.0d0*Datan(1.0d0)
character(len=100) :: par_file
logical :: output_exists

call cpu_time(t1)

alphaf = 0.d0
do ii=1, a2-a1+1
	exponent = ii-1
        alphaf = alphaf+a(part,ii)*order**(exponent)
end do

deltaf = 0.d0
do ii=1, d2-d1+1
	exponent=ii-1
        deltaf = deltaf+d(part,ii)*order**(exponent)
end do

thetaf = 0.d0
do ii=1, th2-th1+1
	exponent = ii-1
        thetaf = thetaf+t(part,ii)*order**(exponent)
end do

gammaf = 0.d0
do ii=1, g2-g1+1
	exponent = ii-1
        gammaf = gammaf+g(part,ii)*order**(exponent)
end do

epsilonf = 0.d0
do ii=1, e2-e1+1
	exponent = ii-1
        epsilonf = epsilonf+eps(part,ii)*order**(exponent)
end do

lambda_cf = x(1)
max_blazef = x(2)

i=cur_dat
j=cur_ord
chi2=0.d0

! point by point
do k=1, n_pnt
	
	!determine the value of the blaze function in the spectrum i, the order j and the point k
	blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, i, j, k)
	if (teluric(i,j,k)) then
		
		result_wav(i,j,k)=observed_wav(i,j,k)
		result_int(i,j,k)=observed_int(i,j,k)/blaze_value
		
		observed_sigma(i,j,k) = 0.01
		
		! determine and add chi2 for this point
		chi2 = chi2+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2

		if (debug_chi2_func_presx) then
			write (201,*) "prsx:", i, j, k, blaze_value, observed_wav(i,j,k), observed_int(i,j,k), synthetic_int(i,j,k),result_int(i,j,k),&
			observed_sigma(i,j,k), alphaf, deltaf, thetaf, gammaf, epsilonf, chi2
		end if
	else
		result_wav(i,j,k) = 0
		result_int(i,j,k) = 0
		observed_sigma(i,j,k) = 0
	end if	
	
end do ! k = n_pnt

! printing current iteration to chi2_func_presx.tmp
if (tmp_chi2_func_presx) then
	if(r==2) then
		call cpu_time(t2)
		write(54, *) x, chi2, t2-t1
		r=0
	else
		r=r+1
	end if
end if

chi2_func_presx = chi2
return

end function chi2_func_presx

end module chi2_func_presx_module
