! teluric.f90
! module for bsimplex
! module creates logical array so that any point inside any teluric line can be ignored

! Pavel Dolezal September 26th 2023
! pavel.dolezal@matfyz.cuni.cz

module teluric_rem_module

contains 

subroutine teluric_rem

use loader_module
use variables_module

implicit none

double precision :: telx1, telx2, tely1, tely2, new_corr, old_corr
character(len=30) :: file_name, corrname
character(len=120) :: line
integer :: ioerr, iocorrerr, i, j, k, n

print *, "Teluric module started..."
n_telpoint = 0

! open file with heliocentric corrections 
inquire (file = "common_input/heliocorr.asc",  exist = file_existence)
if (file_existence) then
	open (20, file = "common_input/heliocorr.asc")
else
	print *, "Error: No common_input/heliocorr.asc file"
	stop
end if

! open file with teluric spectrum
inquire (file = "common_input/teluric.asc",  exist = file_existence)
if (file_existence) then
	open (80,file = "common_input/teluric.asc")
else
	print *, "Error: No common_input/teluric.asc. file."
	stop
end if

print *, "file common_input/teluric.asc opened"

if (debug_teluric) then
        open (456, file = "debug_teluric.tmp")
end if

Do i=1, n_dat
	print *, "Teluric correction of dataset #", i
	write (unit = file_name, FMT = "(A,A,I5.5,A)") "input_rs3/",trim(prefix), list(i), ".asc"
	print *, file_name
	flush(6)
	
	! find correction for current dataset
	do while (.not.(file_name(11:) == corrname))
		50 read (20, "(A)", iostat = ioerr)  line
		if (line (1:1) == "#") go to 50 ! commentary detection
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
			! check if datapoint is in any teluric line
			telx1=0
			do
			200 	read (80, "(A)", iostat = ioerr) line
				if ((line (1:1) == "#").or.(len_trim(line) == 0)) go to 200
				if (ioerr /= 0) then
					print *, "Error EOF common_input/teluric.asc"
					stop
				end if
				
				read (line, *) telx2, tely2
				telx2 = telx2*(new_corr/299792.458+1) ! aply heliocentric correction to point of teluric spectrum
				
				if (telx2 > observed_wav(i,j,k)) then
					if (telx1 < 0.01) then
						print *, "Error: 1st point of common_input/teluric.asc > 1st data point; ",&
						"(i, j, k, telx1) = ", i, j, k, telx1
						stop
					else
						exit
					end if
				end if
				telx1 = telx2
				tely1=tely2
			end do
			
			backspace(80)
			
			if (tely1<teltreshold) then
				teluric(i,j,k) = .false.
				if (debug_teluric) then
					write (456,*) i,j,k,observed_wav(i,j,k)
				end if
				n_telpoint = n_telpoint + 1
			end if
			
		end  do ! k = n_pnt
		
		do n = 1, 1000
			backspace (80)
		end do
		
	end do ! j = n_ord
	
	rewind (80)
	
end do ! i = n_dat

close (80)
close (20)
close (456)

return

end subroutine teluric_rem

end module teluric_rem_module
