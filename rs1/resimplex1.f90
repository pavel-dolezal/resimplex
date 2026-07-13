! resimplex1.f90
!
! Pavel Dolezal September 2nd 2025
! pavel.dolezal@matfyz.cuni.cz
!
! 1st part of the resimplex software suite intended for automatically rectify echelle spectra
! this part takes gamma, epsilon, theta, delta, alpha, lambda_c and max_blaze as constants for each order of each dataset
! i.e. in every order - these constants are completely independent of other orders and datasets
! results are these constants which are meant to be fitted with polynomials in resimplex2
!
!
! the blaze function
! lambda - wavelength
!			 [sin {pi * beta(lambda) * X_r(lambda)}]^2
! R_r(lambda) = A_r = ------------------------------------------------
!			   [pi * beta(lambda) * X_r(lambda)]^2
!	
! X_r(lambda) = r_e * {1-(lambda / lambda_c,r)}
!
! r_e - echelle order
! r - our order system 
! r is shifted echelle order for convenience so that the order with the shortest central wavelength is numbered r = 1 
! r = 126 - r_e
!
! beta_r(lambda) = gamma_r * 10^-10 * (lambda_c,r - epsilon_r - lambda)**4 + theta_r * 10^-8 (lambda_c,r - epsilon_r - lambda)^3
!		    + delta_r * 10^-6 * (lambda_c,r - epsilon_r - lambda)**2 + alpha_r
!
! gamma_r = g_5*r^4 + g4*r^3 + g3*r^2 + g2*r + g1 
! epsilon_r = e_5*r^4 + e4*r^3 + e3*r^2 + e2*r + e1
! theta_r = t_5*r^4 + t4*r^3 + t3*r^2 + t2*r + t1
! delta_r = d_5*r^4 + d4*r^3 + d3*r^2 + d2*r + d1
! alpha_r = a_5*r^4 + a4*r^3 + a3*r^2 + a2*r + a1
!
! lambda_c,r - central wavelength of the order r 
! lambda_c,r = l9*r^8 + ... + l2*r + l1
!
! A_r - local maximum of the blaze function within the order r, i.e. max_blaze
! A_r = m9*r^8 + ... + m2*r + m1
!
! for full documentation see: https://

program resimplex1

use loader_module
use chi2_func_module
use amoeba_mod
use heliocorr_module
use teluric_rem_module
use variables_module

implicit none

double precision :: ftol, e_max, e_lambda(8), e_lc
integer :: ITMAX

! z ... chi2 corresponding to the simplex p
! p ... simplex
! x ... array of parameters
! chi2 ... chi2 for current order
integer :: i, ioerr, iter, j, posmin, it
double precision :: alpha, delta, theta, gamma, epsilon
character(len=300) :: line
double precision :: chi2, ftol_array(2)
double precision, allocatable :: x(:), e(:), z(:), p(:,:), max_blaze(:), lambda_c(:)
character (len=100) :: prefix2

! # of datasets, # of orders, # of points
print *, "n_dat, n_ord, n_pnt:"
4 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 4
read (line,*, err=990, end=990) n_dat, n_ord, n_pnt
print *, n_dat, n_ord, n_pnt

allocate (observed_wav(n_dat, n_ord, n_pnt))
allocate (observed_int(n_dat, n_ord, n_pnt))
allocate (observed_sigma(n_dat, n_ord, n_pnt))
allocate (synthetic_wav(n_dat, n_ord, n_pnt))
allocate (synthetic_int(n_dat, n_ord, n_pnt))
allocate (result_int(n_dat, n_ord, n_pnt))
allocate (result_wav(n_dat, n_ord, n_pnt))
allocate (list(n_dat))
allocate (teluric(n_dat, n_ord, n_pnt))
allocate (lambda_c(n_ord))
allocate (max_blaze(n_ord))

teluric(:,:,:) = .true.

print *, "prefix of data files: "
6 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 6
read (line,*, err=990, end=990) prefix
print *, prefix

print *, "gain of CCD: "
7 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 7
read (line,*, err=990, end=990) gain
print *, gain

print *, "The highest number in data files names: "
8 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 8
read (line,*, err=990, end=990) n_high
print *, n_high

print *, "Final order to compute total chi2: "
9 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 9
read (line,*, err=990, end=990) end_chi2
print *, end_chi2

print *, "rectify?"
50 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 50
read (line,*, err=990, end=990) rectify
print *, rectify

print *, "are synthetic spectra available?"
read (*,*, err=990, end=990) synth_avail
print *, synth_avail

print *, "debug loader?"
read (*,*, err=990, end=990) debug_loader
print *, debug_loader

print *, "debug chi2_func?"
read (*,*, err=990, end=990) debug_chi2_func
print *, debug_chi2_func

print *, "tmp_file?"
read (*,*, err=990, end=990) tmp_file
print *, tmp_file

print *, "debug heliocorr?"
read (*,*, err=990, end=990) debug_heliocorr
print *, debug_heliocorr

print *, "debug teluric?"
read (*,*, err=990, end=990) debug_teluric
print *, debug_teluric

15 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 15

print *, "Initial parameters:"
read(line, *, err=990, end=990) alpha, delta, theta, gamma, epsilon
read(*,*, err=990, end=990) lambda_c
read(*,*, err=990, end=990) max_blaze(1)

max_blaze(:44) = max_blaze(1)
max_blaze(45:) = max_blaze(1)/10

print *, "alpha = ", alpha
print *, "delta = ", delta
print *, "theta = ", theta
print *, "gamma = ", gamma
print *, "epsilon = ", epsilon
print *, "lambda_c = ", lambda_c
print *, "max blaze = ", max_blaze

! number of free parameters
ndim = 5+2*n_dat

print *, "ndim = ", ndim

allocate (x(ndim))
allocate (e(ndim))
allocate (p(ndim+1, ndim))
allocate (z(ndim+1))

print *, "Initial steps e() for alpha, delta, theta, gamma, epsilon:"
55 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 55
read (line,*, err=990, end=990) e(1:5) 
print *, e(1:5)

print *, "Initial steps e() for lambda_c:"
read (*,*, err=990, end=990) e_lc 
e(6:6+n_dat-1) = e_lc 
print *, e(6:6+n_dat-1)

print *, "Initial steps e() for max_blaze:"
read (*,*, err=990, end=990) e_max 
e(6+n_dat:6+2*n_dat-1) = e_max 
print *, e(6+n_dat:6+2*n_dat-1)

print *, "The teluric threshold:"
60 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 60
read (line,*, err=990, end=990) teltreshold
print *, teltreshold

print *, "Array of fractional convergences:"
62 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 62
read (line,*, err=990, end=990) ftol_array
print *, ftol_array

print *, "ITMAX:"
64 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 64
read (line,*, err=990, end=990) ITMAX
print *, ITMAX

print *, "Output comment:"
70 read (*, "(A)", err=990, end=990) comment
if ((comment(1:1) == "#").or.(len_trim(comment) == 0)) goto 70
print *, comment

! load observed and synthetic data
call loader

! remove teluric lines
call teluric_rem

first_run = .true.

open(20, file = "output/input_for_rs2/alpha.tmp")
open(21, file = "output/input_for_rs2/delta.tmp")
open(23, file = "output/input_for_rs2/theta.tmp")
open(26, file = "output/input_for_rs2/gamma.tmp")
open(27, file = "output/input_for_rs2/epsilon.tmp")
open(24, file = "output/input_for_rs2/lambda_c.tmp")
open(25, file = "output/input_for_rs2/max_blaze.tmp")

print *, "Open chi2_func.tmp"
open(50, file = "chi2_func.tmp", status = "unknown", position = "append")

open(51, file = "blaze1.tmp")
write (51, *) "# order & wavelength [A] & blaze_value (max_blaze = constant) & dataset "
open(52, file = "blaze2.tmp")
write (52, *) "# order & wavelength [A] & blaze_value (max_blaze = variable) & dataset "

iter = 0 

! open files to save rectified spectra
if (rectify) then
	do i=1, n_dat	
		n_unit = 300 + i
		write(out_name, "(A,A,I5.5,A)") "output/",trim(prefix), list(i), ".asc"
		print *, out_name, "opened"
		open(unit = n_unit, file = out_name, status = "replace")
		write (n_unit,*) "# wavelength [A] & normalised_intensity I_lambda [] & sigma [] & interpolated? [logical] & blaze_value [] & ",&
		"order_parity []"

		write (n_unit,*) "#",comment
	end do
end if

rectifinal = .false.

print *, "Start iterations"

chi2_all = 0.d0

! file for chi2 per order
if (rectify) then
	open (15,file="chi2_order.tmp")
end if

DO cur_ord=1, n_ord
	print *
	print *, "Current order: ", cur_ord
	
	chi2_ord = 0.d0

	! coloring of orders for highlighting the overlap
	color_index = (-1)**cur_ord
	if (color_index<0) then
		color = 1
	else
		color = 2
	end if

	! initialization of x
	x(1) = alpha
	x(2) = delta
	x(3) = theta
	x(4) = gamma
	x(5) = epsilon
	x(6:6+n_dat-1) = lambda_c(cur_ord)
	x(6+n_dat:6+2*n_dat-1) = dabs(max_blaze(cur_ord))
	
	print *, "x = ", x
	print *, "e = ", e
	print *, "chi2(x) = ", chi2_func(x)
	
	! multi-iteration
	DO it=1,2
		ftol = ftol_array(it)
		do i=1, ndim+1
			p(i,:)=x
			if (i<ndim+1) p(i,i)=p(i,i)+e(i)
		end do
	
		print *, "The initial simplex for it ", it, ":"
		do i=1, ndim+1  
			write(*, *) "p(",i,",:) = ", p(i,:)
		end do

		! determine chi2 of the initial simplex
		do i=1, ndim+1
			x=p(i,:)
			z(i) = chi2_func(x)
			print *, "z(", i, ") = ", z(i)
		end do
				
		r = 0 ! reset counter for chi2_func

		print *, "current order, current iteration: ", cur_ord, it
		flush(6)

		! let's roll
		call amoeba(p,z,ftol,chi2_func,iter,ITMAX)

		! determine the best result
		posmin = minloc(z,1)
		x = p(posmin,:)
		chi2 = z(posmin)
		
		print *, "x = ", x
		print *, "chi2 = ", chi2
	end do ! it = 2

	print *, "Results for order: ", cur_ord
	print *, "x = ", x
	print *, "chi2 = ", chi2

	! total chi2
	if (cur_ord<end_chi2) then
		chi2_all = chi2_all+chi2
	end if

	write (20,*) cur_ord, abs(x(1)) !alpha
	write (21,*) cur_ord, x(2) !delta
	write (23,*) cur_ord, x(3) !theta
	write (26,*) cur_ord, x(4) !gamma
	write (27,*) cur_ord, x(5) !epsilon
	write (24,*) cur_ord, x(6:5+n_dat) !lambda_c
	write (25,*) cur_ord, dabs(x(6+n_dat:5+2*n_dat)) !max_blaze

	! call chi2_func for final rectification, compute chi2 with proper sigma_i and for the whole spectra
	if (rectify) then
		rectifinal = .true.
		chi2 = chi2_func(x)
		rectifinal = .false.
	end if

	if (rectify) then 
		write (15,*) cur_ord, chi2_ord
	end if			

END DO ! cur_ord=n_dat

close(15)

close(20)
close(21)
close(22)
close(23)
close(25)
close(26)
close(27)

close(50)
close(60)

print *, "chi2_all = ", chi2_all
print *, "# of teluric points removed = ", n_telpoint
print *

deallocate (x)
deallocate (e)
deallocate (p)
deallocate (z)
deallocate (observed_wav)
deallocate (observed_int)
deallocate (observed_sigma)
deallocate (synthetic_wav)
deallocate (synthetic_int)
deallocate (result_int)
deallocate (result_wav)
deallocate (list)
deallocate (teluric)
deallocate (max_blaze)
deallocate (lambda_c)

if (rectify) then
	do i=1, n_dat	
		close(unit = 300+i)
	end do
end if

stop

! errors
990 continue
write (*,*) "Error during parameters reading!"

deallocate (x)
deallocate (e)
deallocate (p)
deallocate (z)
deallocate (observed_wav)
deallocate (observed_int)
deallocate (observed_sigma)
deallocate (synthetic_wav)
deallocate (synthetic_int)
deallocate (result_int)
deallocate (result_wav)
deallocate (list)
deallocate (teluric)

if (rectify) then
	do i=1, n_dat
		close(300+i)
	end do
end if

end program resimplex1
