obj = \
variables.o \
heliocorr.o \
loader.o \
blaze.o \
chi2_func.o \
nrtype.o \
nrutil.o \
amoeba.o \
prefit_a_func.o \
prefit_a.o \
prefit_d_func.o \
prefit_d.o \
prefit_e_func.o \
prefit_e.o \
prefit_t_func.o \
prefit_t.o \
prefit_g_func.o \
prefit_g.o \
teluric.o \
chi2_func_presx.o \
presimplex.o \
run.o \

#opt = -O0 -g -fcheck=all -fbacktrace -Wextra
opt = -O2

all: resimplex2

resimplex2: resimplex2.f90 $(obj)
	gfortran $(opt) $(obj) -o $@ $<

$(obj): %.o:%.f90
	gfortran $(opt) -c -o $@ $<

clean:
	rm *.o *.mod

reset:
	make clean
	make all
