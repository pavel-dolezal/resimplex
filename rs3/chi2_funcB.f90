! chi2_funcB.f90
! module for bsimplex program
! computation of chi2 for given parameters x(...) of the blaze function
! where x contains coefficients of max_blaze polynomial
!
! Pavel Dolezal October 6th 2025
! pavel.dolezal@matfyz.cuni.cz

module	chi2_funcB_module

double precision :: blaze_value, alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, sigma, t1, t2

contains

double precision function chi2_funcB(x)

use loader_module
use blaze_module
use variables_module

implicit none
double precision, dimension(:), intent(in) :: x
double precision, parameter:: max_blaze_const = 1e5
double precision :: chi2, exponent
integer :: i, j, k, ii

call cpu_time(t1)

chi2 = 0.d0
chi2_allB = 0.d0

do j=1, n_ord
	chi2_ord=0.d0
	if (((skipped_orders(j) == 0).or.(rectifinal)).and.((j>min2).and.(j<=max2))) then
		order = j ! convert integer -> real
		
	        color_index = (-1)**j
	        if (color_index < 0 ) then
		        color = 1
	        else
		        color = 2
	        end if

		i=dataset

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
		
		lambda_cf = 0.d0
                do ii=1,size(lc,3)
                        exponent = ii-1
                        lambda_cf = lambda_cf+lc(part,dataset,ii)*order**(exponent)
                end do

		max_blazef = 0.d0
		do ii=1, size(x)
			exponent = ii-1
			max_blazef = max_blazef+x(ii)*order**(exponent)
		end do		

		if (rectifinal) then
			write (32,*) i,j,max_blazef
		end if
		
		if (max_blazef < 0.0) then
			chi2 = chi2 + 1000*abs(max_blazef)
		end if

		do k=1, n_pnt	!across points
			
			!determine the value of the blaze function in the spectrum i, the order j and point k
			blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, i, j, k)
			
			if ((teluric(i,j,k)).or.((rectifinal).and.(wtelur))) then
				
				! divide the value of the flux in this point by the value of the blaze function in this point
				result_wav(i,j,k)=observed_wav(i,j,k)
				result_int(i,j,k)=observed_int(i,j,k)/blaze_value
				
				!observed_sigma(i,j,k) = 0.01
				observed_sigma(i,j,k) = sqrt(gain/observed_int(i,j,k))*result_int(i,j,k)		
				
				! rectification
				if (rectifinal) then
					write (300, *) result_wav(i,j,k), result_int(i,j,k), observed_sigma(i,j,k), skipped_orders(j), blaze_value, color
					if (teluric(i,j,k)) then
						chi2_ord = chi2_ord+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2
					end if
				end if 
				
				! determine and add chi2 for this point
				if (teluric(i,j,k)) then
					chi2 = chi2+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2
				end if							
	
				if (debug_chi2_funcB) then
					write (201,*) i, j, k, blaze_value, observed_wav(i,j,k), observed_int(i,j,k), synthetic_int(i,j,k),result_int(i,j,k),&
					observed_sigma(i,j,k), alphaf, deltaf, thetaf, gammaf, epsilonf, chi2
				end if
			else
				result_wav(i,j,k) = 0
				result_int(i,j,k) = 0
				observed_sigma(i,j,k) = 0
			end if
			
			! print value of blaze function with proper max_blaze and with max_blaze = 100 000
			if (rectifinal) then
				write (52, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
				blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blaze_const, lambda_cf, i, j, k)
				write (51, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
			end if
			
		end do ! k = n_pnt

		if ((rectifinal).and.(j<end_chi2)) then
			chi2_ord2(j) = chi2_ord2(j)+chi2_ord
		end if
	end if

        ! chi2_all based on orders 
        if ((rectifinal).and.(j<end_chi2)) then
                chi2_allB = chi2_allB + chi2_ord
        end if

end do ! j=n_ord

! printing every 10th iteration to chi2_funcB.tmp
if (tmp_chi2_funcB) then
	if (r == 10) then	
		call cpu_time(t2)
		write(55, *) x, chi2, t2-t1
		r = 0
	else
		r = r+1
	end if
end if

chi2_funcB = chi2
return

end function chi2_funcB

end module chi2_funcB_module
