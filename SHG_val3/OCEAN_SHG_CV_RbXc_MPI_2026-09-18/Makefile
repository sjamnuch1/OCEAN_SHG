FC = gfortran
MPIFC = mpifort
FFLAGS = -O2 -std=f2008 -Wall -Wextra -fcheck=bounds

all: shg_rbxc_cv.x

# Default production executable: MPI distributed over conduction bands.
shg_rbxc_cv.x: shg_rbxc_cv_mpi.f90
	$(MPIFC) $(FFLAGS) -o $@ $<

# Retain the original serial implementation as a reference and small-test path.
shg_rbxc_cv_serial.x: shg_rbxc_cv.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f shg_rbxc_cv.x shg_rbxc_cv_serial.x
.PHONY: all clean
