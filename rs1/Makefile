obj = \
variables.o \
heliocorr.o \
loader.o \
blaze.o \
chi2_func.o \
nrtype.o \
nrutil.o \
amoeba.o \
teluric.o \

opt = -O3

all: resimplex1

resimplex1: resimplex1.f90 $(obj)
	gfortran $(opt) $(obj) -o $@ $<

$(obj): %.o:%.f90
	gfortran $(opt) -c -o $@ $<

fresh:
	rm *.out *.err *.tmp *.png || true
clean:
	rm *.o *.mod

reset:
	make fresh
	make clean
	make all
