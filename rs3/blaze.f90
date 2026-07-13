! blaze.f90
! module for resimplex program
! computation of blaze function

! Pavel Dolezal October 6th 2025
! pavel.dolezal@matfyz.cuni.cz
module blaze_module

contains

double precision function blaze(alphab, deltab, thetab, gammab, epsilonb, max_blazeb, lambda_cb, i, j, k)

use loader_module
use variables_module

implicit none

double precision, intent(IN) :: alphab, deltab, thetab, gammab, epsilonb, max_blazeb, lambda_cb
integer, intent(in) :: i, j, k
integer :: echelle_order, our_order
double precision :: X_func, beta, argument, inner
double precision, parameter :: pi = 4.0d0*Datan(1.0d0)

!blaze = max_blazeb * (sin(pi*X_func*beta))^2/(pi*X_func*beta)^2

our_order = j 
echelle_order = 126-our_order

X_func= (echelle_order)*(1-(lambda_cb/observed_wav(i,j,k)))

if (X_func==0) then
	blaze = dabs(max_blazeb)
else
	inner = lambda_cb-epsilonb-observed_wav(i,j,k)
	beta = gammab*1.0d-10*inner*inner*inner*inner+thetab*1.0d-8*inner*inner*inner+deltab*1.0d-6*inner*inner+alphab

	if (beta==0) then
		blaze=dabs(max_blazeb)
	else 
		argument = pi*beta*X_func
		blaze = dabs(max_blazeb)*((sin(argument))/(argument))**2
	end if
end if

return

end function blaze

end module blaze_module
