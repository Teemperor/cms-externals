exampledir = $(top_srcdir)/TauSpinner/examples/CP-tests

# No default example programs
# Only available if proper paths are present
example_PROGRAMS = 

INCLUDES   = -I$(prefix)/include
AM_LDFLAGS = -R $(prefix)/lib

LDADD =                          \
$(FLIBS)                         \
$(prefix)/lib/libTauolaTauSpinner.so \
$(prefix)/lib/libTauolaCxxInterface.so  \
$(prefix)/lib/libTauolaFortran.so \
-L$(with_lhapdf)/lib -lLHAPDF

if HAS_HEPMC
if HAS_MCTESTER

# HepMC
AM_LDFLAGS += -R $(HEPMC_DIR)/lib 
INCLUDES   += -I$(HEPMC_DIR)/include
LDADD      += -L$(HEPMC_DIR)/lib -lHepMC

# MC-TESTER
AM_LDFLAGS += -R $(MCTESTER_DIR)/lib
INCLUDES   += -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD      += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)

example_PROGRAMS += CP-test.exe

CP_test_exe_SOURCES = CP-test.cxx \
                      ../../src/read_particles_from_TAUOLA.cxx

if HAS_PYTHIA8

# pythia8
AM_LDFLAGS += -R $(PYTHIA8_DIR)/lib/archive
INCLUDES   += -I$(PYTHIA8_DIR)/include
LDADD      += -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lpythia8tohepmc

example_PROGRAMS += generate-datafiles/generate.exe

generate_datafiles_generate_exe_SOURCES = generate-datafiles/generate.cxx

endif # PYTHIA8

endif # MC-TESTER
endif # HepMC
