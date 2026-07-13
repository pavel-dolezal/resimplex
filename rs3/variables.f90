module variables_module

integer :: n_dat, n_ord, n_pnt !# of spectra, # of orders, # of points
logical :: debug_loader, debug_chi2_funcA, debug_chi2_funcB, debug_chi2_func2
logical :: tmp_chi2_funcA, tmp_chi2_funcB, tmp_chi2_func2, tmp_prefitA_func, tmp_prefitlc_func
logical :: rectify, rectifinal, synth_debug, debug_heliocorr, wtelur, debug_teluric
logical :: synth_avail

integer :: n_skip, cur_dat, ITMAX1, ITMAX2, max_ord, lc_ord, color, color_index
integer :: numline, r, n_unit, n_telpoint, cur_ord, max2, min2, end_chi2, n_high
double precision :: ftol_pres, order, chi2_all, chi2_all2, chi2_all_ls, chi2_all_ms, gain, chi2_diff

double precision, allocatable :: a(:,:), d(:,:), t(:,:), g(:,:), eps(:,:)
double precision :: chi2_allA, chi2_allB, chi2_all1, chi2_ord, conv_treshold, em_presx, elc_presx
logical :: file_existence, file_selector, maxb_asc_existence
logical, allocatable :: teluric(:,:,:)
character (len=100) :: out_name
character (len=300) :: comment
character (len=10) :: prefix
integer, allocatable :: skipped_orders(:)
integer :: dataset, conv, n_conv, part, split(3), adtge_poly
double precision, allocatable :: lf(:), mf(:), maxb(:,:,:), elc(:), emb(:), lc(:,:,:), lc2(:,:), mb(:,:,:), m2(:,:), &
				free_lc(:,:), free_mb(:,:), chi2_ord2(:)

double precision, allocatable :: observed_wav(:,:,:), observed_int(:,:,:), observed_sigma(:,:,:), &
result_wav(:,:,:), result_int(:,:,:), synthetic_wav(:,:,:), synthetic_int(:,:,:)
double precision :: teltreshold
integer, allocatable :: list(:), skip(:)


end module variables_module
