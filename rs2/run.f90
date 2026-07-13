! run.f90
! module for finding lambda_c and max_blaze for each order and each dataset
! lambda_c and max_blaze are free, meaning they are constant for each order and dataset
! and are not bind by any function
!
! Pavel Dolezal February 19th 2025
! pavel.dolezal@matfyz.cuni.cz

module run_module

contains

subroutine run(min,max)

use variables_module
use loader_module
use chi2_func_module
use amoeba_mod
use presimplex_module

implicit none

double precision, allocatable :: x(:), e(:), z(:), p(:,:)
double precision :: chi2, chi2A, chi2B
integer :: i, it, n, iter, posmin, num
integer, intent(IN) :: min, max

iter = 0

min2=min
max2=max

allocate (x(ndim))
allocate (e(ndim))
allocate (p(ndim+1, ndim))
allocate (z(ndim+1))

! initialization of x
x(a1:a2) = a(part,:)
x(d1:d2) = d(part,:)
x(th1:th2) = t(part,:)
x(g1:g2) = g(part,:)
x(e1:e2) = eps(part,:)

e(a1:a2) = ea
e(d1:d2) = ed
e(th1:th2) = et
e(g1:g2) = eg
e(e1:e2) = eeps

print *, "run...", part
print *, "x = ", x
print *, "e = ", e

print *, "chi2(x) = ", chi2_func(x)
flush(6)
it=0
chi2_all1=0.d0
chi2A=0.d0
chi2B=0.d0
chi2_diff=0.d0
num=0
r=0

DO
	! create initial simplex
	
	do i=1, ndim+1
		p(i,:)=x
		if (i<ndim+1) p(i,i)=p(i,i)+e(i)
	end do

	print *, "The initial simplex:"
	do i=1, ndim+1  
		write(*, *) "p(",i,",:) = ", p(i,:)
	end do

	! determine chi2 of the initial simplex
	do i=1, ndim+1
		x=p(i,:)
		z(i) = chi2_func(x)
		print *, "z(", i, ") = ", z(i)
		flush(6)
	end do
	
	r = 0 ! reset counter for chi2_func
	n = 0 ! reset counter of non skipped orders for chi2_func	
	
	print *, "run ",part, ", convergence: ", it+1, ", previous chi2A-chi2B = ", chi2_diff

	flush(6)

	! let's roll
	call amoeba(p,z,ftol,chi2_func,iter,ITMAX)

	! determine the best result
	posmin = minloc(z,1)
	x = p(posmin,:)
	chi2 = z(posmin)
	print *, "current convergence, chi2: ", it+1, chi2
	print *, "x = ", x

	a(part,:) = x(a1:a2)
	d(part,:) = x(d1:d2)
	t(part,:) = x(th1:th2)
	g(part,:) = x(g1:g2)	
	eps(part,:) = x(e1:e2)
	
        ! stop the loop if final convergence has reduced chi2 by less than conv_treshold
        chi2B = chi2
        it=it+1
        if ((chi2A>10).and.((chi2A-chi2B)<conv_treshold)) then
                print *, "the final convergence #", it
                print *, "chi2A-chi2B = ", chi2A-chi2B
                exit
        else
                chi2_diff = chi2A-chi2B
                chi2A = chi2B
        end if

	if (num==5) then
		call presimplex
		num=0
	else
		num = num+1
	end if

end do ! conv_treshold

! call chi2_func for final rectification, compute chi2 with proper sigma_i
if (rectify) then
	rectifinal = .true.
	chi2 = chi2_func(x)
	rectifinal = .false.
end if

print *, "Results, for run ", part, ":"
print *, "x = ", x
print *, "chi2 = ", chi2

print *, "# of teluric points removed = ", n_telpoint
print *
flush(6)

write (55,*) a(part,:)
write (55,*) d(part,:)
write (55,*) t(part,:)
write (55,*) g(part,:)
write (55,*) eps(part,:)

deallocate (x)
deallocate (e)
deallocate (p)
deallocate (z)

end subroutine run

end module run_module
