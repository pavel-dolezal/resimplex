#!/usr/bin/gawk -f

BEGIN{
  c = 299792459.;  # m/s
  Ang = 1.e-10;  # m
#  sigma_Int = 0.01;  # assumed uncertainty

  print "# JD & lambda [m] & I_lambda [] normalised intensity & sigma_I [] & dataset & filename & interpolation [0/1] & value of blaze []"

#  str = "4000-4025 4097-4155 4260-4280 4308-4405 4450-4490 4700-4725 4817-4935 5008-5025 6670-6690";
#  str = "4000-4025 4097-4155 4260-4280 4308-4405 4450-4490 4700-4725 4817-4935 5008-5025 6500-6690";	# Halpha
#  str = "4125-4155 4260-4280 4380-4405 4450-4490 4700-4725 4912-4935 5008-5025 6670-6690";	# without Balmer
#  str = "4125-4155 4260-4280 4380-4405 4450-4490 4700-4725 4912-4935 5008-5025 6500.6-6501.6 6509.9-6511.4 6520.5-6522.3 6525.3-6530.0 6538.2-6541.7 6546.2-6547.2 6549.7-6552.3 6554.7-6557.1 6559.2-6563.8 6565.1-6571.7 6573.2-6574.3 6576.0-6580.5 6588.0-6598.7 6603.5-6690";	# without Balmer, with Halpha, without telluric
#  str = "5008-5025";	# HeI 5015
#  str = "3990-4320 4360-4840 4880-5450";  # no Hbeta, no Hgamma
#  str = "4504-5423";
  str = "4504-8950";
#  str = "4504-5495";
  n = split(str, region, "  *");
  for (i=1; i<=n; i++){
    j = split(region[i], l, "-")
    reg1[i] = l[1];
    reg2[i] = l[2];
    print "# reg[", i, "] = ", reg1[i], " to ", reg2[i], " Ang"
  }
}
(ARGIND==1) && (FNR>=1){
  filename = $1;
  jd[filename] = 2400000.0 + $2;
# sigma_Int__[filename] = 0.01;
#  weight = $3;
#  sigma_Int__[filename] = 0.01*1./weight;
#  vhelio[filename] = $3*1.e3;  # m/s
}
(ARGIND!=1) && !/^#/ && !/^Exported/{
  gsub("\r", "");
  lambda = $1+0.0;
  Int = $2;
  sigma_Int = $3;
  selected = 0;
  idataset = ARGIND-1;
  interpolation = $4;
  blaze_value = $5;
  for (i=1; i<=n; i++) {
    if ((lambda >= reg1[i]) && (lambda <= reg2[i]) && (Int > 0.0)) {
      selected = 1;
    }
  }

  # use only Ondrejov dataset for H_alpha
#  if ((lambda >= 6500) && (lambda <= 6690) && (idataset < 33)) {
#   selected = 0;
#  }

  if (selected == 1) {
    lambda = $1*Ang;
#    Delta_lambda = -lambda * vhelio[filename]/c;  # this is NOT needed!
    Delta_lambda = 0.;
   # sigma_Int = sigma_Int__[FILENAME]
    dataset = FILENAME;

    printf("%.6f  %.8e  %.8f  %.8f  %d  %s  %s  %.8f\n", jd[FILENAME], lambda + Delta_lambda, Int, sigma_Int, idataset, dataset, interpolation, blaze_value);
  }
}

