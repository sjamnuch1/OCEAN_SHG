# Unified build for the OCEAN SHG core-valence sector.
#
# Put this Makefile in the same directory as:
#   shg_gmres.f90
#   shg_valence_source.f90
#   shg_xc_minus_valence.f90
#   shg_rbxc_cv_mpi.f90
#   shg_core_source.f90
#   shg_cv_manifest.f90
#   shg_y1_cv.f90
#   shg_t1_gv_contract.f90
#
# By default that directory is assumed to be directly below the OCEAN source
# root.  Override OCEAN_ROOT when it is elsewhere, for example:
#
#   make -j8 OCEAN_ROOT=/path/to/OCEAN
#
# OCEAN itself must be compiled first so that its .o and .mod files exist.

OCEAN_ROOT ?= ..
ARCH_MAKEFILE ?= $(OCEAN_ROOT)/Makefile.arch

# Import the compiler, optimization, FFT, BLAS, and OpenMP settings used to
# build OCEAN.  The leading '-' lets `make help` and the standalone targets be
# inspected even before OCEAN_ROOT has been set correctly.
-include $(ARCH_MAKEFILE)

MPIFORT ?= ftn
OPTIONS ?= -O2
LDFLAGS ?=
FFTWL ?=
BLAS ?=
OPENMP ?=
EXTRA_FFLAGS ?=
EXTRA_LDFLAGS ?=

OCEAN_DIR := $(OCEAN_ROOT)/OCEAN2
MODULES_DIR := $(OCEAN_ROOT)/Modules
SCREEN_DIR := $(OCEAN_ROOT)/SCREEN/src
CNBSE_DIR := $(OCEAN_ROOT)/CNBSE

OCEAN_FFLAGS := $(OPTIONS) \
	-I$(OCEAN_DIR) \
	-I$(MODULES_DIR) \
	-I$(SCREEN_DIR) \
	$(EXTRA_FFLAGS)

STANDALONE_FFLAGS := $(OPTIONS) $(EXTRA_FFLAGS)

# This is the same OCEAN object set used by the individual core-valence
# Makefiles.  OCEAN.o is intentionally absent because it defines ocean.x's
# main program and would conflict with the SHG program entry points.
OCEAN_OBJECTS := \
	$(OCEAN_DIR)/AI_kinds.o \
	$(OCEAN_DIR)/OCEAN_mpi.o \
	$(OCEAN_DIR)/OCEAN_system.o \
	$(OCEAN_DIR)/OCEAN_bloch.o \
	$(OCEAN_DIR)/OCEAN_obf.o \
	$(OCEAN_DIR)/OCEAN_multiplet.o \
	$(OCEAN_DIR)/long_range.o \
	$(OCEAN_DIR)/OCEAN_psi.o \
	$(OCEAN_DIR)/OCEAN_energies.o \
	$(OCEAN_DIR)/getabb.o \
	$(OCEAN_DIR)/getomega.o \
	$(OCEAN_DIR)/gamfcn.o \
	$(OCEAN_DIR)/jimel.o \
	$(OCEAN_DIR)/nbsemkcmel.o \
	$(OCEAN_DIR)/intval.o \
	$(OCEAN_DIR)/newgetylm.o \
	$(OCEAN_DIR)/newgetprefs.o \
	$(OCEAN_DIR)/newthreey.o \
	$(OCEAN_DIR)/elsdch.o \
	$(OCEAN_DIR)/cinv.o \
	$(OCEAN_DIR)/OCEAN_timekeeper.o \
	$(OCEAN_DIR)/OCEAN_invdrv.o \
	$(OCEAN_DIR)/OCEAN_val_energy.o \
	$(OCEAN_DIR)/OCEAN_read_tmels.o \
	$(OCEAN_DIR)/OCEAN_val_states.o \
	$(OCEAN_DIR)/OCEAN_bubble.o \
	$(OCEAN_DIR)/OCEAN_get_rho.o \
	$(OCEAN_DIR)/optim.o \
	$(OCEAN_DIR)/OCEAN_hyb_louie_levine.o \
	$(OCEAN_DIR)/OCEAN_ladder.o \
	$(OCEAN_DIR)/OCEAN_rixs_holder.o \
	$(OCEAN_DIR)/OCEAN_corewidths.o \
	$(OCEAN_DIR)/OCEAN_action.o \
	$(OCEAN_DIR)/OCEAN_filenames.o \
	$(OCEAN_DIR)/OCEAN_WRR.o \
	$(OCEAN_DIR)/OCEAN_fxc.o \
	$(OCEAN_DIR)/OCEAN_load_data.o \
	$(MODULES_DIR)/FFT_wrapper.o \
	$(MODULES_DIR)/jlmfft.o \
	$(CNBSE_DIR)/zcain/zobj/getprefs.o \
	$(CNBSE_DIR)/zcain/zobj/getylm.o \
	$(CNBSE_DIR)/zcain/zobj/fancyvector.o \
	$(SCREEN_DIR)/ocean_sphericalharmonics.o

XC_OBJECTS := \
	shg_gmres.o \
	shg_valence_source.o \
	shg_xc_minus_valence.o

Y1_OBJECTS := \
	shg_gmres.o \
	shg_core_source.o \
	shg_cv_manifest.o \
	shg_y1_cv.o

PROGRAMS := \
	shg_xc_minus_valence.x \
	shg_rbxc_cv_mpi.x \
	shg_y1_cv.x \
	shg_t1_gv_contract.x

SOURCES := \
	shg_gmres.f90 \
	shg_valence_source.f90 \
	shg_xc_minus_valence.f90 \
	shg_rbxc_cv_mpi.f90 \
	shg_core_source.f90 \
	shg_cv_manifest.f90 \
	shg_y1_cv.f90 \
	shg_t1_gv_contract.f90

.DEFAULT_GOAL := all

.PHONY: all ocean standalone serial-rbxc check-sources check-ocean \
	print-config help clean

all: check-sources $(PROGRAMS)

# Programs that link to OCEAN's BSE machinery.
ocean: shg_xc_minus_valence.x shg_y1_cv.x

# Programs whose source has no OCEAN module dependency.
standalone: shg_rbxc_cv.x shg_t1_gv_contract.x

check-sources:
	@missing=0; \
	for source in $(SOURCES); do \
		if test ! -f "$$source"; then \
			echo "ERROR: missing source file $$source"; \
			missing=1; \
		fi; \
	done; \
	test $$missing -eq 0

check-ocean:
	@test -f "$(ARCH_MAKEFILE)" || { \
		echo "ERROR: cannot find $(ARCH_MAKEFILE)"; \
		echo "Set OCEAN_ROOT to the compiled OCEAN source root."; \
		exit 1; \
	}
	@test -f "$(OCEAN_DIR)/ocean_system.mod" || { \
		echo "ERROR: $(OCEAN_DIR)/ocean_system.mod is missing."; \
		echo "Build OCEAN in OCEAN_ROOT before building the GMRES programs."; \
		exit 1; \
	}
	@test -f "$(MODULES_DIR)/FFT_wrapper.o" || { \
		echo "ERROR: $(MODULES_DIR)/FFT_wrapper.o is missing."; \
		echo "Build the OCEAN support modules first."; \
		exit 1; \
	}

# Shared OCEAN-based modules.
shg_gmres.o: shg_gmres.f90 | check-ocean
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

shg_valence_source.o: shg_valence_source.f90 | check-ocean
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

shg_core_source.o: shg_core_source.f90 | check-ocean
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

shg_cv_manifest.o: shg_cv_manifest.f90 | check-ocean
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

# X_c^-(omega) from the valence source.
shg_xc_minus_valence.o: shg_xc_minus_valence.f90 \
		shg_gmres.o shg_valence_source.o
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

shg_xc_minus_valence.x: $(XC_OBJECTS) | check-ocean
	$(MPIFORT) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ \
		$(XC_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

# MPI R_cv application, distributed over conduction bands.
shg_rbxc_cv_mpi.o: shg_rbxc_cv_mpi.f90
	$(MPIFORT) $(STANDALONE_FFLAGS) -c -o $@ $<

shg_rbxc_cv_mpi.x: shg_rbxc_cv_mpi.o
	$(MPIFORT) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ $^ \
		$(OPENMP) $(OPTIONS)

# Y_bc(2 omega) core-BSE solve using the dumped R_b X_c vectors.
shg_y1_cv.o: shg_y1_cv.f90 \
		shg_gmres.o shg_core_source.o shg_cv_manifest.o
	$(MPIFORT) $(OCEAN_FFLAGS) -c -o $@ $<

shg_y1_cv.x: $(Y1_OBJECTS) | check-ocean
	$(MPIFORT) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ \
		$(Y1_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

# Final d_a^\dagger Y_bc contraction; this source is standalone and serial,
# but MPIFORT keeps the compiler/runtime consistent with the other programs.
shg_t1_gv_contract.o: shg_t1_gv_contract.f90
	$(MPIFORT) $(STANDALONE_FFLAGS) -c -o $@ $<

shg_t1_gv_contract.x: shg_t1_gv_contract.o
	$(MPIFORT) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ $^ \
		$(OPENMP) $(OPTIONS)

# Optional build of the older serial R_cv implementation.  It is not part of
# `all`; invoke `make serial-rbxc` when shg_rbxc_cv.f90 is present.
serial-rbxc: shg_rbxc_cv_serial.x

shg_rbxc_cv_serial.o: shg_rbxc_cv.f90
	$(MPIFORT) $(STANDALONE_FFLAGS) -c -o $@ $<

shg_rbxc_cv_serial.x: shg_rbxc_cv_serial.o
	$(MPIFORT) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ $^ \
		$(OPENMP) $(OPTIONS)

print-config:
	@echo "OCEAN_ROOT       = $(OCEAN_ROOT)"
	@echo "ARCH_MAKEFILE    = $(ARCH_MAKEFILE)"
	@echo "MPIFORT          = $(MPIFORT)"
	@echo "OPTIONS          = $(OPTIONS)"
	@echo "EXTRA_FFLAGS     = $(EXTRA_FFLAGS)"
	@echo "EXTRA_LDFLAGS    = $(EXTRA_LDFLAGS)"

help:
	@echo "Core-valence SHG build targets"
	@echo ""
	@echo "  make [-jN] [OCEAN_ROOT=/path/to/OCEAN]  Build all four programs"
	@echo "  make ocean                              Build the two OCEAN-linked solvers"
	@echo "  make standalone                         Build RbXc and final contraction"
	@echo "  make serial-rbxc                        Build the older serial RbXc code"
	@echo "  make print-config                       Show the selected build settings"
	@echo "  make clean                              Remove generated local build files"

clean:
	$(RM) \
		shg_gmres.o \
		shg_valence_source.o \
		shg_xc_minus_valence.o \
		shg_rbxc_cv_mpi.o \
		shg_rbxc_cv_serial.o \
		shg_core_source.o \
		shg_cv_manifest.o \
		shg_y1_cv.o \
		shg_t1_gv_contract.o \
		shg_gmres.mod \
		shg_valence_source.mod \
		shg_core_source.mod \
		shg_rbxc_cv_manifest_io.mod \
		$(PROGRAMS) \
		shg_rbxc_cv_serial.x

