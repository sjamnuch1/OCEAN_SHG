include ../Makefile.arch

OCEAN_DIR = ../OCEAN2
FLAGS = $(OPTIONS) -I$(OCEAN_DIR) -I../Modules -I../SCREEN/src

OCEAN_OBJECTS = \
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
	../Modules/FFT_wrapper.o \
	../Modules/jlmfft.o \
	../CNBSE/zcain/zobj/getprefs.o \
	../CNBSE/zcain/zobj/getylm.o \
	../CNBSE/zcain/zobj/fancyvector.o \
	../SCREEN/src/ocean_sphericalharmonics.o

Y1_CV_OBJECTS = shg_gmres.o shg_core_source.o shg_cv_manifest.o shg_y1_cv.o

.PHONY: all ocean_objects clean

all: shg_y1_cv.x shg_t1_gv_contract.x

ocean_objects:
	@test -f $(OCEAN_DIR)/ocean.x || { \
		echo "Build the unmodified OCEAN source from its repository root first."; \
		exit 1; \
	}
	@test -f ../Modules/FFT_wrapper.o || { \
		echo "Missing OCEAN support objects; run make in the repository root."; \
		exit 1; \
	}

shg_gmres.o: shg_gmres.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_core_source.o: shg_core_source.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_cv_manifest.o: shg_cv_manifest.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_y1_cv.o: shg_y1_cv.f90 shg_gmres.o shg_core_source.o shg_cv_manifest.o
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_y1_cv.x: $(Y1_CV_OBJECTS) ocean_objects
	$(MPIFORT) $(LDFLAGS) -o $@ $(Y1_CV_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

shg_t1_gv_contract.x: shg_t1_gv_contract.f90
	$(MPIFORT) $(OPTIONS) -o $@ $<

clean:
	$(RM) $(Y1_CV_OBJECTS) shg_y1_cv.x shg_t1_gv_contract.x \
		shg_gmres.mod shg_core_source.mod shg_rbxc_cv_manifest_io.mod
