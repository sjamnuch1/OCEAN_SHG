include ../Makefile.arch

OCEAN_DIR = ../OCEAN2
MODULE_DIR = ../Modules

FLAGS = $(OPTIONS) -I$(MODULE_DIR) -I$(OCEAN_DIR) -I../SCREEN/src

# Build against the object set from the original OCEAN Makefile.  The link
# deliberately omits OCEAN.o because it contains ocean.x's program entry point.
# No SHG-specific source or object is required under OCEAN2/.
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

XC_OBJECTS = shg_gmres.o shg_core_source.o shg_xc_minus.o
RBXC_OBJECTS = shg_core_vertex.o shg_rbxc_core.o
RBXC_REALSPACE_OBJECTS = shg_core_vertex.o shg_core_vertex_realspace.o \
	shg_rbxc_core_realspace.o

.PHONY: all ocean_objects clean

all: shg_xc_minus.x shg_rbxc_core.x shg_rbxc_core_realspace.x

ocean_objects:
	@test -f $(OCEAN_DIR)/ocean.x || { \
		echo "Build the unmodified OCEAN source from its repository root first."; \
		exit 1; \
	}
	@test -f ../Modules/FFT_wrapper.o || { \
		echo "Missing original OCEAN support objects; run 'make' in the repository root."; \
		exit 1; \
	}

shg_gmres.o: shg_gmres.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_core_source.o: shg_core_source.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_xc_minus.o: shg_xc_minus.f90 shg_gmres.o shg_core_source.o
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_xc_minus.x: $(XC_OBJECTS) ocean_objects
	$(MPIFORT) $(LDFLAGS) -o $@ $(XC_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

shg_core_vertex.o: shg_core_vertex.f90 ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_rbxc_core.o: shg_rbxc_core.f90 shg_core_vertex.o ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_rbxc_core.x: $(RBXC_OBJECTS) ocean_objects
	$(MPIFORT) $(LDFLAGS) -o $@ $(RBXC_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

shg_core_vertex_realspace.o: shg_core_vertex_realspace.f90 shg_core_vertex.o ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_rbxc_core_realspace.o: shg_rbxc_core_realspace.f90 \
		shg_core_vertex.o shg_core_vertex_realspace.o ocean_objects
	$(MPIFORT) $(FLAGS) -c -o $@ $<

shg_rbxc_core_realspace.x: $(RBXC_REALSPACE_OBJECTS) ocean_objects
	$(MPIFORT) $(LDFLAGS) -o $@ $(RBXC_REALSPACE_OBJECTS) $(OCEAN_OBJECTS) \
		$(FFTWL) $(BLAS) $(OPENMP) $(FFTWL) $(OPTIONS)

clean:
	$(RM) $(XC_OBJECTS) $(RBXC_OBJECTS) $(RBXC_REALSPACE_OBJECTS) \
		shg_gmres.mod shg_core_source.mod shg_core_vertex.mod \
		shg_core_vertex_realspace.mod shg_xc_minus.x shg_rbxc_core.x \
		shg_rbxc_core_realspace.x
