! heliocorr.f90
! module for resimplex 
! module aplies new heliocentric correction to observed and synthetic data

! Pavel Dolezal September 26th 2023
! pavel.dolezal@matfyz.cuni.cz

module heliocorr_module

contains

subroutine heliocorr

use variables_module

implicit none

character(len=30) :: file_name, corrname
character(len=120) :: line
integer :: ioerr, iocorrerr, i, j, k

double precision :: old_obs, old_synth, new_corr, old_corr ! [km/s]
double precision, parameter :: c = 299792.458 ! [km/s]

! open file with heliocentric corrections 
inquire (file = "common_input/heliocorr.asc",  exist = file_existence)
if (file_existence) then
	open (20, file = "common_input/heliocorr.asc")
else
	print *, "Error: No common_input/heliocorr.asc file"
	stop
end if

do i = 1, n_dat
	print *, "Heliocentric correction of dataset #", i
	write (unit = file_name, FMT = "(A,A,I5.5,A)") "input_rs3/", trim(prefix), list(i), ".asc"

	! find correction for current dataset
	do while (.not.(file_name(11:) == corrname))
		50 read (20, "(A)", iostat = ioerr)  line
		if (line (1:1) == "#") go to 50
		if (ioerr /= 0) then
			print *, "Error: EOF in common_input/heliocorr.asc"
			stop
		end if
		read (line, *, iostat = iocorrerr) corrname, old_corr, new_corr
		if (iocorrerr /= 0) then
			print *, "Error during correction read, i = ", i
			stop
		end if
	end do
	
	do j=1, n_ord
		do k=1, n_pnt
			! apply correction to observed data
			! old correction -> laboratory value -> new correction
			if (debug_heliocorr) then
                                old_obs = observed_wav(i,j,k)
                        end if

                        observed_wav(i,j,k) = observed_wav(i,j,k)/(old_corr/c+1)
                        observed_wav(i,j,k) = observed_wav(i,j,k)*(new_corr/c+1)

                        if (debug_heliocorr) then
                                write (345,*) i,j,k,"obs: ", old_obs, " -> ", observed_wav(i,j,k)
                        end if
		end do ! k
	end do	! j
end do	! i

close (20)

print *, "helliocorr ends succesfully"
flush(6)

return

end subroutine heliocorr

end module heliocorr_module
