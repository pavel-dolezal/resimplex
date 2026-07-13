! chi2_func.f90
! module for resimplex program
! computation of chi2 for given parameters x(...) of the blaze function
! x contains coefficients of polynomial parameters alpha, delta, epsilon, theta, gamma
!
! Pavel Dolezal September 10th 2025
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
double precision :: chi2, exponent, chi2_ord

logical :: skipping
double precision, parameter :: pi = 4.0d0*Datan(1.0d0)
integer :: i, j, k, n_unit2

call cpu_time(t1)

! chi2 computation
chi2 = 0.0d0

! order by order
do j=1, n_ord
	chi2_ord = 0.d0
	if (((skipped_orders(j) == 0).or.(rectifinal)).and.((j>min2).and.(j<=max2))) then
		order = j ! convert integer -> real

                color_index = (-1)**j
                if (color_index < 0 ) then
                        color = 1
                else
                        color = 2
                end if
		
		alphaf = 0.d0
		do i=a1, a2
			exponent = i-a1
			alphaf = alphaf+x(i)*order**(exponent)
		end do
		
		deltaf = 0.d0
		do i=d1, d2
			exponent = i-d1
			deltaf = deltaf+x(i)*order**(exponent)
		end do
		
		thetaf = 0.d0
		do i=th1, th2
			exponent = i-th1
			thetaf = thetaf+x(i)*order**(exponent)
		end do

		gammaf = 0.d0
		do i=g1, g2
			exponent = i-g1
			gammaf = gammaf+x(i)*order**(exponent)
		end do

		epsilonf = 0.d0
		do i=e1, e2
			exponent = i-e1
			epsilonf = epsilonf+x(i)*order**(exponent)
		end do
		
		! dataset by dataset
		do i=1, n_dat
			lambda_cf = lc(i,j)
			max_blazef = m(i,j)
			
			! during final rectification write parameters to variables.tmp
			if ((i==1).and.(rectifinal)) then
				write (15, *) j, alphaf, deltaf, thetaf, gammaf, epsilonf, lambda_cf, max_blazef
			end if

			! point by point
			do k=1, n_pnt
				
				!determine the value of the blaze function in the spectrum i, the order j and the point k
				blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blazef, lambda_cf, i, j, k)
				
				if ((teluric(i,j,k)).or.((rectifinal).and.(wtelur))) then
					
					result_wav(i,j,k)=observed_wav(i,j,k)
					result_int(i,j,k)=observed_int(i,j,k)/blaze_value
					
					observed_sigma(i,j,k) = sqrt(gain/observed_int(i,j,k))*result_int(i,j,k)		

					! determine and add chi2 for this point
					chi2 = chi2+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2
					chi2_ord = chi2_ord+((result_int(i,j,k)-synthetic_int(i,j,k))/observed_sigma(i,j,k))**2

					! rectification
					if (rectifinal) then
						n_unit2 = 300+i
						write (n_unit2, *) result_wav(i,j,k), result_int(i,j,k), observed_sigma(i,j,k), skipped_orders(j), blaze_value, color
					end if 

					if (debug_chi2_func) then
						write (200,*) i, j, k, blaze_value, observed_wav(i,j,k), observed_int(i,j,k), synthetic_int(i,j,k),result_int(i,j,k),&
						observed_sigma(i,j,k), alphaf, deltaf, thetaf, gammaf, epsilonf, chi2
					end if
				else
					result_wav(i,j,k) = 0
					result_int(i,j,k) = 0
					observed_sigma(i,j,k) = 0
				end if	
				
				! print value of blaze function with max_blaze = 100 000
				if (rectifinal) then
					write (52, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
					blaze_value = blaze(alphaf, deltaf, thetaf, gammaf, epsilonf, max_blaze_const+0*j, lambda_cf, i, j, k)
					write (51, *) j, observed_wav(i,j,k), observed_wav(i,j,k)-lambda_cf, blaze_value, i
				end if
				
			end do ! k = n_pnt
		end do ! i = n_dat
	        ! chi2_all1 based on orders without Paschen series
	        if (rectifinal) then
			write (56,*) j, chi2_ord
			if (j<end_chi2) then
				chi2_all1 = chi2_all1 + chi2_ord
			end if		
		end if
	end if
end do ! j=n_ord

if (chi2<0.001) then
	print *, "Error(chi2_func): chi2 = ", chi2
	stop
end if

! print current iteration to chi2_func.tmp
if (tmp_chi2_func) then
	if (r==10) then
		call cpu_time(t2)
		write(50, *) x, chi2, t2-t1
		r=0
	else
		r=r+1
	end if
end if

chi2_func = chi2
return

end function chi2_func

end module chi2_func_module
