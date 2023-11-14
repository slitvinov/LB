.POSIX:
.SUFFIXES:
.SUFFIXES: .f

FC = gfortran
FCFLAGS = -O2 -g

all: bgk2
.f:
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -f bgk2
bgk2: bgk2.par
