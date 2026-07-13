module variables_module

integer :: n_dat, n_ord, n_pnt !# of spectra, # of orders, # of points
logical :: debug_loader, debug_teluric, rectify, rectifinal, first_run, tmp_file, debug_chi2_func, varblaze, debug_heliocorr
logical :: synth_avail

integer :: ndim, cur_dat, end_chi2
integer :: run, numline, r, n_unit, n_telpoint, cur_ord, n, color, color_index, n_high
double precision :: chi2_ord, chi2_all, gain
logical :: variable(5), file_existence, file_selector
logical, allocatable :: teluric(:,:,:)
character (len=100) :: out_name
character (len=300) :: comment
character (len=10) :: prefix

double precision, allocatable :: observed_wav(:,:,:), observed_int(:,:,:), observed_sigma(:,:,:), &
result_wav(:,:,:), result_int(:,:,:), synthetic_wav(:,:,:), synthetic_int(:,:,:), blaze_par(:,:,:)
double precision :: teltreshold
integer, allocatable :: list(:)

end module variables_module
