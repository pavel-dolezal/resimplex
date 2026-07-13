! chi2_func.f90
! module for resimplex program
! computation of chi2 for given parameters x(...) of the blaze function

! Pavel Dolezal September 4th 2025
! pavel.dolezal@matfyz.cuni.cz
module	chi2_func_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function chi2_func(x)

use loader_module
use blaze_module
use variables_module

implicit none
double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, l(8), order

double precision, parameter :: pi = 4.0d0*Datan(1.0d0)
integer :: i, j, k, n_unit2

call cpu_time(t1)

! chi2 computation
chi2 = 0.0d0

! determine chi2
j=cur_ord

do i=1, n_dat
	
	alphaf = x(1)
	deltaf = x(2)
	thetaf = x(3)
	gammaf = x(4)
	epsilonf = x(5)	

	lambda_cf = x(5+i)
	max_blazef = dabs(x(5+n_dat+i))
	
	if (max_blazef == 0.0) then
		max_blazef = 1.0
	end if
		
	do k=1, n_pnt	!across points
		
		!determine the value of the blaze function in the spectrum i, the order j and point k
		blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, i, j, k)
		
		if (teluric(i,j,k)) then
			
			! divide the value of the flux in this point by the value of the blaze function in this point
			result_wav(i,j,k)=observed_wav(i,j,k)
			result_int(i,j,k)=observed_int(i,j,k)/blaze_value
			
			observed_sigma(i,j,k) = sqrt(gain/observed_int(i,j,k))*result_int(i,j,k)		
					
			! rectification
			if (rectifinal) then
				chi2_ord = chi2_ord+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2

				n_unit2 = 300+i
				write (n_unit2, *) result_wav(i,j,k), result_int(i,j,k), observed_sigma(i,j,k), "0", blaze_value, color
			end if 
			
			! determine and add chi2 for this point
			chi2 = chi2+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2
			
			if (debug_chi2_func) then
				open (200,file = "chi2_func_debug.tmp", status = "unknown", position = "append")
				write (200,*) i, j, k, blaze_value, observed_wav(i,j,k), observed_int(i,j,k), synthetic_int(i,j,k),result_int(i,j,k),&
				observed_sigma(i,j,k), chi2
				close (200)
			end if
		else
			result_wav(i,j,k) = 0
			result_int(i,j,k) = 0
			observed_sigma(i,j,k) = 0
		end if	
		
		! print value of blaze function with max_blaze = 100 000
		!if ((rectifinal).and.(i==1)) then
		if (rectifinal) then
			write (52, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
			blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blaze_const+0*j, lambda_cf, i, j, k)
			write (51, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
		end if
		
	end do ! k = n_pnt 
end do ! j = n_ord

close(40)

if (rectifinal) then
	first_run = .false.
end if

! printing current iteration to chi2_func.tmp, only every 100th iteration
if (tmp_file) then
	if (r == 100) then	
		call cpu_time(t2)
		write(50, *) x, chi2, t2-t1

		r = 0
	else
		r = r+1
	end if
end if

chi2_func = chi2
return

end function chi2_func

end module chi2_func_module
