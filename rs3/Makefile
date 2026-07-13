obj = \
variables.o \
heliocorr.o \
loader.o \
blaze.o \
chi2_funcA.o \
chi2_funcB.o \
nrtype.o \
nrutil.o \
amoeba.o \
teluric.o  \
prefitlc_func.o \
prefitlc.o \
prefitA_func.o \
prefitA.o \
lambdac_search.o \
maxblaze_search.o \
chi2_func2.o \
presimplex.o \

#opt = -O0 -g -fcheck=all -fbacktrace -Wall -Wextra
opt = -O3

all: resimplex3

resimplex3: resimplex3.f90 $(obj)
	gfortran $(opt) $(obj) -o $@ $<

$(obj): %.o:%.f90
	gfortran $(opt) -c -o $@ $<

fresh:
	rm *.out *.err *.tmp *.o *.mod || true
clean:
	rm *.o *.mod

reset:
	make fresh
	make all
