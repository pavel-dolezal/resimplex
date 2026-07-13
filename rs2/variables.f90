module variables_module

integer :: n_dat, n_ord, n_pnt !# of spectra, # of orders, # of points
logical :: debug_loader, rectify, rectifinal, debug_chi2_func, varblaze, synth_debug, debug_heliocorr, wtelur, &
		debug_chi2_func_presx, tmp_prefit_func, tmp_chi2_func_presx, tmp_chi2_func, debug_teluric, &
		synth_avail
		
integer :: n_skip, n_max, cur_dat, ITMAX, ITMAX1, split(3), adtge_poly, end_chi2, n_high, ndim
integer :: count = 0
integer :: numline, r, n_unit, n_telpoint, cur_ord, part, dataset, fit_min, fit_max
double precision :: chi2_pseudo, ftol_pres, order, chi2_all, conv_treshold, gain, chi2_diff
double precision :: chi2_part(4)
double precision, allocatable :: a(:,:), d(:,:), t(:,:), g(:,:), eps(:,:)
logical :: variable(5), file_existence, file_selector
logical, allocatable :: teluric(:,:,:)
character (len=100) :: out_name
character (len=300) :: comment
character (len=10) :: prefix
integer, allocatable :: skipped_orders(:)
integer :: l1, l2, m1, m2, a1, a2, d1, d2,th1, th2, g1, g2, e1, e2, color, color_index, min2, max2
double precision, allocatable :: lf(:), mf(:), lc(:,:), m(:,:), free_a(:), free_d(:), free_t(:), free_g(:), free_eps(:)
double precision, allocatable :: ea(:), ed(:), et(:), eg(:), eeps(:)
double precision :: ftol, e_max, e_lambda(8), elc, em, chi2_all1, chi2_all2

double precision, allocatable :: observed_wav(:,:,:), observed_int(:,:,:), observed_sigma(:,:,:), &
result_wav(:,:,:), result_int(:,:,:), synthetic_wav(:,:,:), synthetic_int(:,:,:)
double precision :: teltreshold
integer, allocatable :: list(:), skip(:)


end module variables_module
