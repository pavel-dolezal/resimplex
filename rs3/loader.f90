! loader.f90
! module for bsimplex program
! module for loading CHIRON spectra and synthetic spectra into arrays
! expecting n_dat files with spectra each one with n_ord orders
! individual orders in data spectra must be divided by a blank line
! two blank lines at the EOF results in error  "k<size(data,3)"

! Pavel Dolezal June 22nd 2023
! pavel.dolezal@matfyz.cuni.cz

module loader_module

contains

subroutine loader

use heliocorr_module
use variables_module

implicit none

double precision :: jd, lambda, intensity, tol ! variables for synthetic

character (len=30) :: file_name
character (len=100) :: line
integer :: i, j, k, ioerr, dataioerr, n
logical :: file_existence_syn

! maximum difference between synthetic wavelength and observed wavelength [Angstrom]
tol = 1e-3

! for save feature in b_chi2_func.f90
file_selector=.true.

if (debug_loader) then
	open (789,file = "debug_loader.tmp")
end if 

! data loading
i=0
print *, "i, file: "
do n=1, 350
	! load every existing data file in "data" directory
	write (unit = file_name, FMT = "(A,A,I5.5,A)") "input_rs3/", trim(prefix), n, ".txt"
	inquire (file = file_name,  exist = file_existence)
	if (file_existence) then
		open (30, file = file_name)
		i=i+1
		
		!write down used data file number for later rectification
		list(i) = n 
		
		! is # of spectra > # of expected spectra
		if (i > size(observed_wav,1)) then	
			print *, "Error: i = ", i, " > ", "size(observed_wav,1) = ", size(observed_wav,1)
			stop
		end if
		
		print *, i, file_name
		j=0
		ioerr = 0
		
		do ! load current data file order by order
			if (ioerr /= 0) then
				if (j < size(observed_wav,2)) then
					print *, " i = ", i, ", j_max = ", j, "< size(observed_wav,2) = ", size(observed_wav,2)	! # of orders in data file is smaller then expected
					stop
				end if
				exit
			end if
				
			j=j+1
			k=1
			if (debug_loader) then
				write (789,*) "i, j, k, observed_wav(i,j,k), observed_int(i,j,k): "
			end if
			
			do ! point by point
			100 read (30, "(A)", iostat = ioerr)  line
				if (line (1:1) == "#") go to 100 ! commentary or header detection
				
				if ((len_trim(line) == 0).or.(ioerr /= 0)) then	! detection of EOF or an empty line = end of an order
					if ((k-1 < size(observed_wav,3)).and.(k-1 /= 0)) then	! # of points is smaller then expected
						print *, "i = ", i, ", j = ", j, ", k_max = ", k-1, " < SIZE(observed_wav,3) = ", size(observed_wav,3)
						stop
					end if
					exit
				end if
				
				if (j > size(observed_wav,2)) then ! check # of orders
					print *, "Error: i = ", i, ", j = ", j, " > ", "size(observed_wav,2) = ", size(observed_wav,2)
					stop
				end if
				
				if (k > size(observed_wav,3)) then ! check # of points
					print *, "Error: i = ", i, ", j = ", j, ", k = ", k, " > ", "size(observed_wav,3) = ", size(observed_wav,3)
					stop
				end if
				
				! actual load
				read (line, *, iostat=dataioerr) observed_wav(i, j, k), observed_int(i, j, k)
				if (dataioerr /= 0) then
					print *, "Error during data read: i, j ,k = ", i, j, k
					stop
				end if
			
                                if (observed_int(i,j,k) < 0.0) then
                                        print *, "Error: observed intensity < 0; i,j,k: ", i,j,k
                                        stop
                                end if

				! print loaded data
				if (debug_loader) then
					write (789,*) i, j, k, observed_wav(i, j, k), observed_int(i, j, k) 
				end if
				
				k=k+1
			end do
		end do
		close (30)
	end if
end do

if (i < size(observed_wav,1)) then ! # of spectra si smaller then expected
	print *, "i_max = ", i, "< size(observed_wav,1) = ", size(observed_wav,1)
	stop
end if
print *, "# of loaded spectra: ", i

! apply corrected (by teluric lines) heliocentric correction to observed spectra
! synthetic spectra are already corrected!
call heliocorr

! synthetic data loading
! synthetic data are not divided into orders!
do i = 1, size(observed_wav,1)	!across all spectra
	if (synth_avail) then
		inquire (file = "common_input/synthetic.dat", exist = file_existence)
		if (.not.file_existence) then
			print *, "Error: file common_input/synthetic.dat missing"
			stop
		end if
		
		numline =0

		open (10, file = "common_input/synthetic.dat")
		
		write (unit = file_name, FMT = "(A,A,I5.5,A)") "input_rs3/", trim(prefix), list(i), ".asc"	
		print *, file_name
		
		do	j = 1, size(observed_wav,2)	!across all orders
			if (debug_loader) then
				write (789, *) "i, j, k, synthetic_wav(i,j,k,1), synthetic_int(i,j,k,2): "
			end if
			
			do k=1, size(observed_wav,3) !across all data points

				if (observed_wav(i,j,k) < 0.01) then	! check if data point is in a teluric line
					synthetic_wav(i,j,k) = 0
					synthetic_int(i,j,k) = 0
					if (debug_loader) then
						write (789,*) i, j, k, synthetic_wav(i,j,k), synthetic_int(i,j,k)
					end if
				else
				
				!scan the synthetic data for the right point
				do	
				70	read (10, "(A)", iostat = ioerr)  line
					numline = numline+1
					if (line (1:1) == "#") go to 70
					if (ioerr /= 0) then
						print *, "Error: no synth point for i, j, k, observed_wav = ", i, j, k, observed_wav(i,j,k)
						stop
					end if
					
					read (line, *,iostat=dataioerr) JD, lambda, intensity, dataset
					!if ((dataset /= i).or.((len_trim(line)) == 0))go to 70
					if (len_trim(line) == 0) go to 70
					
					lambda = lambda*10d0**10d0
					
					if (dataioerr /= 0) then
						print *, "Error: common_input/synthetic.dat, line = ", numline
						stop
					end if
					
					! check if the synthetic_wav(i, j, k) would match the observed_wav(i, j, k)
					if (abs(lambda-observed_wav(i,j,k))<=tol) then
						synthetic_wav(i, j, k) = lambda
						synthetic_int(i, j, k) = intensity
						
						if (debug_loader) then
							write (789,*) i, j, k, synthetic_wav(i,j,k), synthetic_int(i,j,k)
						end if
						
						exit
					end if
				end do
				end if
			end do
			
			!rewind synthetic.asc because of orders overlay
			do n=1, 500
				backspace(10)
				numline = numline-1
			end do
		end do	
        else
                do j = 1, size(observed_wav,2)
                        do k = 1, size(observed_wav,3) !across all data points
                                synthetic_wav(i,j,k) = observed_wav(i,j,k)
                                synthetic_int(i,j,k) = 1.0
                        end do
                end do
        end if
end do

close (10)

end subroutine loader

end module loader_module
