! resimplex3.f90
!
! 3rd and final part of the resimplex software suite intended for automatic rectification of echelle spectra
!
! This part applies alpha, delta, epsilon, theta, gamma as polynomial functions from resimplex2
! and computes lambda_c, and max_blaze as polynomial functions for each dataset individually. 
!
! z ... chi2 corresponding to the simplex p
! p ... simplex
! x ... array of parameters
!
! Pavel Dolezal October 9th 2025
! pavel.dolezal@matfyz.cuni.cz

program resimplex3

use loader_module
use chi2_funcA_module
use chi2_funcB_module
use amoeba_mod
use heliocorr_module
use teluric_rem_module
use variables_module
use lambdac_search_module
use maxblaze_search_module
use presimplex_module

implicit none

double precision :: ftol

integer :: i, ioerr, iter, j, k, posmin, ord, it, jj, n, ndim
integer, allocatable :: check(:)
character(len=500) :: line
double precision :: chi2_tmp, chi2, chi2_1, chi2_2, chi2_3, chi2_4
double precision, allocatable :: x(:), e(:), z(:), p(:,:)


! # of datasets, # of orders, # of points
print *, "n_dat, n_ord, n_pnt:"
4 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 4
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

print *, "number greater than the highest in data file names:"
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
35 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 35
read (line,*, err=990, end=990) split
print *, split

print *, "adtge poly: "
40 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 40
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

print *, "debug chi2_funcA?"
read (*,*, err=990, end=990) debug_chi2_funcA
print *, debug_chi2_funcA

print *, "debug chi2_funcB?"
read (*,*, err=990, end=990) debug_chi2_funcB
print *, debug_chi2_funcB

print *, "debug chi2_func2?"
read (*,*, err=990, end=990) debug_chi2_func2
print *, debug_chi2_func2

print *, "tmp file for chi2_funcA?"
read (*,*, err=990, end=990) tmp_chi2_funcA
print *, tmp_chi2_funcA

print *, "tmp file for chi2_funcB?"
read (*,*, err=990, end=990) tmp_chi2_funcB
print *, tmp_chi2_funcB

print *, "tmp file for chi2_func2?"
read (*,*, err=990, end=990) tmp_chi2_func2
print *, tmp_chi2_func2

print *, "tmp file for prefitA_func?"
read (*,*, err=990, end=990) tmp_prefitA_func
print *, tmp_prefitA_func

print *, "tmp file for prefitlc_func?"
read (*,*, err=990, end=990) tmp_prefitlc_func
print *, tmp_prefitlc_func

print *, "debug heliocorr?"
read (*,*, err=990, end=990) debug_heliocorr
print *, debug_heliocorr

print *, "order of max_blaze and lambda_c polynomials (max_ord, lc_ord):"
55 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 55
read (line,*, err=990, end=990) max_ord, lc_ord
print *, max_ord, lc_ord

allocate (chi2_ord2(n_ord))
chi2_ord2 = 0.d0

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
allocate (lc(4,n_dat,lc_ord+1))	! coefficients
allocate (mb(4,n_dat,max_ord+1))	! coefficients
allocate (free_mb(n_dat, n_ord))   ! constants from rs2 - to be fitted in fitmb
allocate (free_lc(n_dat,n_ord)) ! constants from rs2 - to be fitted in prefitlc

allocate (elc(lc_ord+1))
allocate (emb(max_ord+1))

allocate (a(4, adtge_poly+1))
allocate (d(4, adtge_poly+1))
allocate (t(4, adtge_poly+1))
allocate (g(4, adtge_poly+1))
allocate (eps(4, adtge_poly+1))

teluric(:,:,:) = .true.
lc = 0.d0
mb = 0.d0
lc(:,:,1) = 4000
mb(:,:,1) = 50000

print *, "Reading coefficients for alpha, delta, theta, gamma, epsilon"
print *, "from input_parameters_for_rs3.txt"

inquire (file = "input_parameters_for_rs3.txt",  exist = file_existence)
if (file_existence) then
        open (25, file = "input_parameters_for_rs3.txt")
else
        print *, "Error: No input_parameters_for_rs3.txt file, run rs2 first"
        stop
end if

60 read (25,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 60
backspace (25)
do i=1,4
	read(25,*, err=990, end=990) a(i,:)
	read(25,*, err=990, end=990) d(i,:)
	read(25,*, err=990, end=990) t(i,:)
	read(25,*, err=990, end=990) g(i,:)
	read(25,*, err=990, end=990) eps(i,:)
end do

close (25)

print *, "Constants:"
print *, "a = ", a
print *, "d = ", d
print *, "t = ", t
print *, "g = ", g
print *, "eps = ", eps

print *, "max_blaze and lambda_c for freeAlc module"
74 read (*,"(A)", err=990, end=990) line
if ((line(1:1)=="#").or.(len_trim(line)==0)) goto 74
read(line, *, err=990, end=990) free_mb(1,1)
free_mb(:,:44) = free_mb(1,1)
free_mb(:,45:) = free_mb(1,1)/10
read (*, *, err=990, end=990) free_lc(1,:)
do i=1,n_dat
	free_lc(i,:) = free_lc(1,:)
end do
print *, "max_blaze: ", free_mb
print *, "lambda_c: ", free_lc

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

! create array of orders with 1 to skip, 0 to not skip
j=1
do i=1, n_ord
	if (i==skip(j)) then
		skipped_orders(i) = 1
		j=j+1
	else
		skipped_orders(i) = 0
	end if
end do !i=n_ord
print *, "array of skipped orders: ", skipped_orders

print *, "Read initial steps e() for lc and m polynomials"
85 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 85
read (line,*, err=990, end=990) elc(:)
read (*,*) emb(:)
read (*,*) elc_presx
read (*,*) em_presx

print *, "elc = ", elc
print *, "emb = ", emb

print *, "The teluric threshold:"
90 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 90
read (line,*, err=990, end=990) teltreshold
print *, teltreshold

print *, "convergence treshold:"
70 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 70
read (line,*, err=990, end=990) conv_treshold
print *, conv_treshold

print *, "# of convergences"
95 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 95
read (line,*, err=990, end=990) n_conv
print *, n_conv

print *, "max # of iterations for presimplex:"
105 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 105
read (line,*, err=990, end=990) ITMAX1
print *, ITMAX1

print *, "max # of iterations for resimplex:"
110 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 110
read (line,*, err=990, end=990) ITMAX2
print *, ITMAX2

print *, "fractional convergence tolerances (ftol):"
112 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 112
read (line,*, err=990, end=990) ftol
print *, ftol

print *, "fractional convergence tolerances for presimplex and prefitA / prefitlc (ftol_pres):"
113 read (*,"(A)", err=990, end=990) line
if ((line(1:1) == "#").or.(len_trim(line)==0)) goto 113
read (line,*, err=990, end=990) ftol_pres
print *, ftol_pres

print *, "output comment"
115 read (*, "(A)", err=990, end=990) comment
if ((comment(1:1) == "#").or.(len_trim(comment) == 0)) goto 115
print *, comment

if (debug_chi2_funcA) then
        open (200, file = "debug_chi2_funcA.tmp")
end if

if (debug_chi2_funcB) then
        open (201, file = "debug_chi2_funcB.tmp")
end if

if (debug_chi2_func2) then
        open (202, file = "debug_chi2_func2.tmp")
end if

if (debug_heliocorr) then
	open (345, file = "debug_heliocorr.tmp")
end if

! load observed and synthetic data
! contains heliocentric correction for observed data
call loader

! remove teluric lines
call teluric_rem

open(51, file = "blaze1.tmp")
write (51, *) "# order & wavelength [A] & blaze_value (max_blaze = constant) & dataset "
open(52, file = "blaze2.tmp")
write (52, *) "# order & wavelength [A] & blaze_value (max_blaze = variable) & dataset "

iter = 0 

rectifinal = .false.


chi2_all = 0.d0
chi2_all1 = 0.d0
chi2_all2 = 0.d0

open (unit=15, file = "chi2_order.tmp")

open (31, file = "output/lc_values_by_order.tmp")
write (31,*) "dataset [] & order [] & lc [A]"
open (32, file = "output/A_values_by_order.tmp")
write (32,*) "dataset [] & order [] & A [ADU]"

open (666, file = "prefit_checklc.tmp")
open (667, file = "prefit_checkA.tmp")

if (tmp_chi2_func2) then
	open (53, file = "chi2_func2.tmp")
end if

if (tmp_chi2_funcA) then
	open (54, file = "chi2_funcA.tmp")
end if

if (tmp_chi2_funcB) then
	open (55, file = "chi2_funcB.tmp")
end if

if (tmp_prefitA_func) then
	open (56, file = "chi2_funcprefitA.tmp")
end if

if (tmp_prefitlc_func) then
	open (57, file = "chi2_funcprefitlc.tmp")
end if


call presimplex

! let's roll
part=1
call lambdac_search(0,split(1))
chi2_all1 = chi2_all1+chi2_all_ls
call maxblaze_search(0,split(1))
chi2_all2 = chi2_all2+chi2_all_ms
chi2_1 = chi2_all2

part=2
call lambdac_search(split(1),split(2))
chi2_all1 = chi2_all1+chi2_all_ls
call maxblaze_search(split(1),split(2))
chi2_all2 = chi2_all2+chi2_all_ms
chi2_2 = chi2_all2-chi2_1

part=3
call lambdac_search(split(2),split(3))
chi2_all1 = chi2_all1+chi2_all_ls
call maxblaze_search(split(2),split(3))
chi2_all2 = chi2_all2+chi2_all_ms
chi2_3 = chi2_all2-chi2_1-chi2_2

part=4
call lambdac_search(split(3),63)
chi2_all1 = chi2_all1+chi2_all_ls
call maxblaze_search(split(3),63)
chi2_all2 = chi2_all2+chi2_all_ms
chi2_4 = chi2_all2-chi2_1-chi2_2-chi2_3

do i=1, n_ord
	write (15,*) i, chi2_ord2(i)
end do

close(15)

print *, "chi2_all_part1 = ", chi2_1
print *, "chi2_all_part2 = ", chi2_2
print *, "chi2_all_part3 = ", chi2_3
print *, "chi2_all_part4 = ", chi2_4	

print *, "chi2_all lambdac_search: ",chi2_all1
print *, "chi2_all maxblaze_search: ", chi2_all2
print *, "# of teluric points removed = ", n_telpoint
print *
flush(6)

open (30, file = "output/poly_coeff_lc.tmp")
write (30,*) "# part [] & dataset [] & poly coeff of lambda_c(r) r^0 & r^1 & r^2 & ... []"
do j=1, 4
	do i=1, n_dat
		write (30,*) j, i, lc(j,i,:) 
	end do
end do
close(30)

open (34, file = "output/poly_coeff_A.tmp")
write (34,*) "# part [] & dataset [] & poly coeff of max_blaze(r) r^0 & r^1 & r^2 & ... [] "
do j=1,4
	do i=1, n_dat
		write (34,*) j, i, mb(j,i,:)
	end do
end do
close(34)

close(31)
close(32)


close(20)
close(21)
close(22)
close(23)
close(25)

close(666)
close(667)

close(200)
close(201)
close(202)

close(53)
close(54)
close(55)
close(56)
close(57)

close(345)

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

close(20)
close(21)
close(22)
close(23)
close(25)

close(53)
close(54)
close(55)
close(56)
close(57)

close(666)
close(667)

close(200)
close(201)
close(202)

close(345)

end program resimplex3
