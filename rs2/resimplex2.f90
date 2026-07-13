! resimplex2.f90
!
! 2nd part of resimplex software suite inteded for automatic rectification of echelle spectra.
!
! This part models gamma, epsilon, theta, delta and alpha as polynomial functions 
! of the variable r (order), these functions are computed once and applied to all datasets.
!
! Lambda_c and A are treated as constants for each order and each dataset, independet of other orders.
! They are taken from the results of the resimplex1 program and periodicaly updated to 
! improve consistency with other parameters during the convergence process.
!
! x ... array of parameters
! p ... simplex
! z ... chi2 corresponding to the simplex p
!
! Pavel Dolezal September 5th 2025
! pavel.dolezal@matfyz.cuni.cz

program resimplex2

use loader_module
use chi2_func_module
use amoeba_mod
use heliocorr_module
use teluric_rem_module
use variables_module
use run_module
use prefit_a_module
use prefit_d_module
use prefit_e_module
use prefit_t_module
use prefit_g_module

implicit none

integer :: i, ioerr, j
integer, allocatable :: check(:)
character(len=500) :: line

! # of datasets, # of orders, # of points
print *, "n_dat, n_ord, n_pnt:"
5 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 5
read (line,*, err=990, end=990) n_dat, n_ord, n_pnt
print *, n_dat, n_ord, n_pnt

print *, "prefix of data files: "
10 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 10
read (line,*, err=990, end=990) prefix
print *, prefix

print *, "gain of CCD: "
15 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 15
read (line,*, err=990, end=990) gain
print *, gain

print *, "The highest number in data files names: "
20 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 20
read (line,*, err=990, end=990) n_high
print *, n_high

print *, "Final order to compute total chi2: "
25 read (*, "(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 25
read (line,*, err=990, end=990) end_chi2
print *, end_chi2

print *, "split: "
30 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 30
read (line,*, err=990, end=990) split
print *, split

print *, "adtge poly: "
35 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 35
read (line,*, err=990, end=990) adtge_poly
print *, adtge_poly

print *, "rectify?"
50 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 50
read (line,*, err=990, end=990) rectify
print *, rectify

print *, "are synthetic spectra available?"
read (*,*, err=990, end=990) synth_avail
print *, synth_avail

print *, "rectify with teluric lines?"
read (*,*, err=990, end=990) wtelur
print *, wtelur

print *, "debug loader?"
read (*,*, err=990, end=990) debug_loader
print *, debug_loader

print *, "debug chi2_func_presx?"
read (*,*, err=990, end=990) debug_chi2_func_presx
print *, debug_chi2_func_presx

print *, "tmp file for chi2_func_presx?"
read (*,*, err=990, end=990) tmp_chi2_func_presx
print *, tmp_chi2_func_presx

print *, "debug chi2_func?"
read (*,*, err=990, end=990) debug_chi2_func
print *, debug_chi2_func

print *, "tmp file for chi2_func?"
read (*,*, err=990, end=990) tmp_chi2_func
print *, tmp_chi2_func

print *, "tmp file for prefit_func?"
read (*,*, err=990, end=990) tmp_prefit_func
print *, tmp_chi2_func

print *, "debug heliocorr?"
read (*,*, err=990, end=990) debug_heliocorr
print *, debug_heliocorr

print *, "debug teluric?"
read (*,*, err=990, end=990) debug_teluric
print *, debug_teluric

allocate (observed_wav(n_dat, n_ord, n_pnt))
allocate (observed_int(n_dat, n_ord, n_pnt))
allocate (observed_sigma(n_dat, n_ord, n_pnt))
allocate (synthetic_wav(n_dat, n_ord, n_pnt))
allocate (synthetic_int(n_dat, n_ord, n_pnt))
allocate (result_int(n_dat, n_ord, n_pnt))
allocate (result_wav(n_dat, n_ord, n_pnt))
allocate (list(n_dat))
allocate (teluric(n_dat, n_ord, n_pnt))
allocate (skipped_orders(n_ord))
allocate (lc(n_dat,n_ord))
allocate (m(n_dat,n_ord))
allocate (free_a(n_ord))
allocate (free_d(n_ord))
allocate (free_t(n_ord))
allocate (free_g(n_ord))
allocate (free_eps(n_ord))

allocate (a(4, adtge_poly+1))
allocate (d(4, adtge_poly+1))
allocate (t(4, adtge_poly+1))
allocate (g(4, adtge_poly+1))
allocate (eps(4, adtge_poly+1))

allocate (ea(size(a,2)))
allocate (ed(size(d,2)))
allocate (et(size(t,2)))
allocate (eg(size(g,2)))
allocate (eeps(size(eps,2)))

teluric(:,:,:) = .true.

print *, "Reading starting coefficients for alpha, delta, theta, gamma, epsilon"
45 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 45
read(line, *, err=990, end=990) a(1,:)
read(*,*, err=990, end=990) d(1,:)
read(*,*, err=990, end=990) t(1,:)
read(*,*, err=990, end=990) g(1,:)
read(*,*, err=990, end=990) eps(1,:)

! all parts of spectra starts with the same polynomial parameters
do i=2,4
	a(i,:) = a(1,:)
	d(i,:) = d(1,:)
	t(i,:) = t(1,:)
	g(i,:) = g(1,:)
	eps(i,:) = eps(1,:)
end do

inquire (file="input_for_rs2/alpha.tmp", exist = file_existence)
if (file_existence) then
	open (14, file="input_for_rs2/alpha.tmp")
else
	print *, "Error: no alpha.tmp file (result of rs1)"
	stop
end if
print *, "Reading alpha.tmp"
do j=1, n_ord
        read (14,*, err=990, end=990) i, free_a(j)
	if (i/=j) then
		print *, "Error: alpha.tmp, i=/j", i, j
		stop
	end if
end do
close (14)

inquire (file="input_for_rs2/delta.tmp", exist = file_existence)
if (file_existence) then
	open (16, file="input_for_rs2/delta.tmp")
else
	print *, "Error: no delta.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/delta.tmp"
do j=1, n_ord
        read (16,*, err=990, end=990) i, free_d(j)
	if (i/=j) then
		print *, "Error: delta.tmp, i=/j", i, j
		stop
	end if
end do
close (16)

inquire (file="input_for_rs2/theta.tmp", exist = file_existence)
if (file_existence) then
	open (17, file="input_for_rs2/theta.tmp")
else
	print *, "Error: no theta.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/theta.tmp"
do j=1, n_ord
        read (17,*, err=990, end=990) i, free_t(j)
	if (i/=j) then
		print *, "Error: theta.tmp, i=/j", i, j
		stop
	end if
end do
close (17)

inquire (file="input_for_rs2/epsilon.tmp", exist = file_existence)
if (file_existence) then
	open (18, file="input_for_rs2/epsilon.tmp")
else
	print *, "Error: no epsilon.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/epsilon.tmp"
do j=1, n_ord
        read (18,*, err=990, end=990) i, free_eps(j)
	if (i/=j) then
		print *, "Error: epsilon.tmp, i=/j", i, j
		stop
	end if
end do
close (18)

inquire (file="input_for_rs2/gamma.tmp", exist = file_existence)
if (file_existence) then
	open (19, file="input_for_rs2/gamma.tmp")
else
	print *, "Error: no gamma.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/gamma.tmp"
do j=1, n_ord
        read (19,*, err=990, end=990) i, free_g(j)
	if (i/=j) then
		print *, "Error: gamma.tmp, i=/j", i, j
		stop
	end if
end do
close (19)

print *, "Reading of lambda_c"

inquire (file="input_for_rs2/lambda_c.tmp", exist = file_existence)
if (file_existence) then
	open (12, file="input_for_rs2/lambda_c.tmp")
else
	print *, "Error: no lambda_c.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/lambda_c.tmp"
do j=1, n_ord
        read (12,*, err=990, end=990) i, lc(:,j)
	if (i/=j) then
		print *, "Error: lambda_c.tmp, i=/j", i, j
		stop
	end if
end do
close (12)

print *, "Reading of max_blaze.tmp"
inquire (file="input_for_rs2/max_blaze.tmp", exist=file_existence)
if (file_existence) then
	open (13, file="input_for_rs2/max_blaze.tmp")
else
	print *, "Error: no max_blaze.tmp file (result of rs1)"
	stop
end if
print*, "Reading input_for_rs2/max_blaze.tmp"
do j=1, n_ord
        read (13,*, err=990, end=990) i, m(:,j)
	if (i/=j) then
		print *, "Error: max_blaze.tmp, i/=j, i,j= ", i,j	
		stop
	end if
end do
close(13)

print *, "Initial parameters:"
print *, "a = ", a
print *, "d = ", d
print *, "t = ", t
print *, "g = ", g
print *, "lc = ", lc
print *, "m = ", m

n_skip = 1
ioerr = 0
print *, "skipped orders:"
57 read (*, "(A)", ERR=990, END=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 57
do while (ioerr==0)
        allocate (check(n_skip))
        read(line,*,iostat=ioerr) check
        if ((ioerr/=0).and.(ioerr/=-1)) then
                print *, "Error skiped orders, line = ", line
                stop
        end if
        n_skip = n_skip+1
        deallocate (check)
end do
n_skip = n_skip-2

print *, "n_skip = ", n_skip
if (n_skip<=0) then
        print *, "Error: no orders to skip, write '99' instead of nothing"
        stop
end if

allocate (skip(n_skip))

read(line, *, err=990, end=990) skip(:)
print *, skip

do i=1, n_skip-1
        if (skip(i)>=skip(i+1)) then
                print *, "Error: skipped orders are not sorted"
                stop
        end if
end do

! create array of orders with 1 to skip, 0 to not skip corresponding order
j=1
do i=1, n_ord
	if (j>size(skip)) then
		exit
	end if
	if (i==skip(j)) then
		skipped_orders(i) = 1
		j=j+1
	else
		skipped_orders(i) = 0
	end if
end do
print *, "skipped orders array (1-skip, 0-do not skip): ", skipped_orders

print *, "Read initial steps e() for alpha, delta, theta, lambda_c, max_blaze"
60 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 60
read (line,*, err=990, end=990) ea
print *, "ea = ", ea
read (*,*) ed
print *, "ed = ",ed
read (*,*) et
print *, "et = ", et
read (*,*) eg
print *, "eg = ",eg
read (*,*) eeps
print *, "eeps = ", eeps
read (*,*) elc
print *, "elc = ", elc
read (*,*) em
print *, "em = ", em

print *, "teluric threshold:"
65 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 65
read (line,*, err=990, end=990) teltreshold
print *, teltreshold

print *, "convergence treshold:"
70 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 70
read (line,*, err=990, end=990) conv_treshold
print *, conv_treshold

print *, "fractional convergence:"
75 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 75
read (line,*, err=990, end=990) ftol
print *, ftol

print *, "maximum of convergences:"
80 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 80
read (line,*, err=990, end=990) ITMAX1
print *, ITMAX1

print *, "maximum of convergences (presimplex):"
85 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 85
read (line,*, err=990, end=990) ITMAX
print *, ITMAX

print *, "rectification comment:"
90 read (*, "(A)", err=990, end=990) comment
if ((comment(1:1) == "#").or.(len_trim(comment) == 0)) goto 90
print *, comment

! load observed and synthetic data
call loader

! remove teluric lines
call teluric_rem

print *, "Open chi2_func.tmp"
open(50, file = "chi2_func.tmp")
open(54, file = "chi2_func_presx.tmp")
open (56, file = "chi2_order.tmp")

open(57, file = "chi2_func_fit.tmp", status = "unknown", position = "append")

open(51, file = "blaze1.tmp")
write (51, *) "# order [] & wavelength [A] & blaze_value [ADU] (max_blaze = constant) & dataset []"
open(52, file = "blaze2.tmp")
write (52, *) "# order [] & wavelength [A] & blaze_value [ADU] (max_blaze = variable) & dataset []"

open (55, file="output/input_parameters_for_rs3.txt")

rectifinal = .false.

! POSITIONS of parameters alpha, delta, theta, epsilon, gamma inside array x
a1=1
a2=size(a,2)
d1=a2+1
d2=a2+size(d,2)
th1=d2+1
th2=d2+size(t,2)
g1=th2+1
g2=th2+size(g,2)
e1=g2+1
e2=g2+size(eps,2)

print *, "adtge 1/2: ", a1,a2,d1,d2,th1,th2,g1,g2,e1,e2

do j=1, n_dat
	if (rectify) then
		n_unit = 300 + j
		write(out_name, "(A,A,I5.5,A)") "output/", trim(prefix), list(j), ".asc"
		print *, out_name, "opened"
		open(unit = n_unit, file = out_name, status = "replace")
		write (n_unit,*) "# wavelength [A] & normalised_intensity I_lambda [] & sigma [] & interpolated? [logical]", &
		" & blaze_value [] & order_parity []"
		write (n_unit,*) "#",comment
	end if
end do

open (unit=15, file = "par_values_by_order.tmp")
write (15,*) "order & alpha & delta & theta & gamma & epsilon & lambda_c & max_blaze (1st dataset only)"

ndim = size(a,2)+size(d,2)+size(t,2)+size(g,2)+size(eps,2)
print *, "ndim = ", ndim

if (debug_chi2_func) then
	open (200,file = "chi2_func_debug.tmp")
	write (200,*) "i & j & k & blaze value & observed wavelength & observed intensity ", &
	"& synthetic intensity & result intensity & sigma & chi2"
end if

if (debug_chi2_func_presx) then
	open (201,file = "chi2_func_presx_debug.tmp")
	write (201,*) "i & j & k & blaze value & observed wavelength & observed intensity ", &
	"& synthetic intensity & result intensity & sigma & chi2"
end if

open (666, file = "prefit_check_a.tmp")
open (667, file = "prefit_check_d.tmp")
open (668, file = "prefit_check_e.tmp")
open (669, file = "prefit_check_t.tmp")
open (670, file = "prefit_check_g.tmp")

open (135, file = "fit_results.tmp")

call prefit_a
call prefit_d
call prefit_e
call prefit_t
call prefit_g

flush(135)
close (135)

chi2_all=0.d0
part = 1
call run(0,split(1))
chi2_all=chi2_all+chi2_all1
chi2_part(1) = chi2_all1

part = 2
call run(split(1),split(2))
chi2_all=chi2_all+chi2_all1
chi2_part(2) = chi2_all1

part = 3
call run(split(2),split(3))
chi2_all=chi2_all+chi2_all1
chi2_part(3) = chi2_all1

part = 4
call run(split(3),n_ord+1)
chi2_all=chi2_all+chi2_all1
chi2_part(4) = chi2_all1

print *, "chi2_part(1) = ", chi2_part(1)
print *, "chi2_part(2) = ", chi2_part(2)
print *, "chi2_part(3) = ", chi2_part(3)
print *, "chi2_part(4) = ", chi2_part(4)
print *, "chi2_all = ", chi2_all

open (70, file="output/rs2_free_max_blaze.tmp")
open (71, file="output/rs2_free_lambda_c.tmp")
do j=1,n_ord
	write (70,*) j, m(:,j)
	write (71,*) j, lc(:,j)
end do

close(15)
close(50)
close(54)
close(55)
close(56)
close(57)
close(60)
close(70)
close(71)
close(200)
close(201)
close(666)
close(667)
close(668)
close(669)
close(670)

do i=1, n_dat
	n_unit = 300+i
	close(n_unit)
end do

deallocate (observed_wav)
deallocate (observed_int)
deallocate (observed_sigma)
deallocate (synthetic_wav)
deallocate (synthetic_int)
deallocate (result_int)
deallocate (result_wav)
deallocate (list)
deallocate (teluric)
!deallocate (max_blaze)
!deallocate (lambda_c)

if (rectify) then
	do i=1, n_dat	
		close(unit = 300+i)
	end do
end if

stop

! errors
990 continue
write (*,*) "Error during parameters reading!"

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

close(15)
close(50)
close(54)
close(55)
close(56)
close(57)
close(60)
close(200)
close(201)
close(666)
close(667)
close(668)
close(669)
close(670)

end program resimplex2
