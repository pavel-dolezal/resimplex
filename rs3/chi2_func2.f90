! chi2_func2.f90
! module for resimplex program
! computation of chi2 for subroutine presimplex and for given parameters x(...) of the blaze function
! x contains lambda_c and max_blaze for given order and dataset
!
! Pavel Dolezal October 6th 2025
! pavel.dolezal@matfyz.cuni.cz

module	chi2_func2_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function chi2_func2(x)

use loader_module
use blaze_module
use variables_module

implicit none

double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, exponent
integer :: i, j, k, ii

call cpu_time(t1)

alphaf = 0.d0
do ii=1, size(a,2)
        exponent = ii-1
        alphaf = alphaf+a(part,ii)*order**(exponent)
end do

deltaf = 0.d0
do ii=1, size(d,2)
        exponent = ii-1
        deltaf = deltaf+d(part,ii)*order**(exponent)
end do

thetaf = 0.d0
do ii=1, size(t,2)
        exponent = ii-1
        thetaf = thetaf+t(part,ii)*order**(exponent)
end do

gammaf = 0.d0
do ii=1, size(g,2)
        exponent = ii-1
        gammaf = gammaf+g(part,ii)*order**(exponent)
end do

epsilonf = 0.d0
do ii=1, size(eps,2)
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

		if (debug_chi2_func2) then
			write (202,*) "prsx:", i, j, k, blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, &
			observed_wav(i,j,k), observed_int(i,j,k), synthetic_int(i,j,k),result_int(i,j,k),&
			observed_sigma(i,j,k), chi2
		end if
	else
		result_wav(i,j,k) = 0
		result_int(i,j,k) = 0
		observed_sigma(i,j,k) = 0
	end if	
	
end do ! k = n_pnt

! printing every 10th iteration to chi2_func2.tmp
if (tmp_chi2_func2) then
	if(r==10) then
		call cpu_time(t2)
		write(53, *) x, chi2, t2-t1
		r=0
	else
		r=r+1
	end if
end if

chi2_func2 = chi2
return

end function chi2_func2

end module chi2_func2_module
