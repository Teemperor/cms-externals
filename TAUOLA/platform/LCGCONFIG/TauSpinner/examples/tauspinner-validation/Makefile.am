exampledir = $(top_srcdir)/TauSpinner/examples/tauspinner-validation

example_PROGRAMS = 

INCLUDES = -I$(prefix)/include        \
           -I$(prefix)/include/Tauola \
           -I$(top_srcdir)/TauSpinner/include/

AM_LDFLAGS = -R $(prefix)/lib -R $(with_lhapdf)/lib

LDADD =                          \
$(FLIBS)                         \
$(prefix)/lib/libTauolaTauSpinner.so \
$(prefix)/lib/libTauolaCxxInterface.so \
$(prefix)/lib/libTauolaFortran.so      \
-L$(with_lhapdf)/lib -lLHAPDF

if HAS_HEPMC

AM_LDFLAGS += -R $(HEPMC_DIR)/lib

INCLUDES += -I$(HEPMC_DIR)/include
LDADD    += -L$(HEPMC_DIR)/lib -lHepMC

if HAS_MCTESTER

AM_LDFLAGS  += -R $(MCTESTER_DIR)/lib

example_PROGRAMS +=  hepmc-tauola-redecay.exe \
                     tauspinner-validation-comparison.exe \
                     tauspinner-validation-plots.exe \
                     test-bornAFB/bornAFB-test.exe \
                     test-ipol/test-ipol.exe

hepmc_tauola_redecay_exe_SOURCES = hepmc-tauola-redecay.cxx \
                                   ../../src/read_particles_from_TAUOLA.cxx

tauspinner_validation_comparison_exe_SOURCES = tauspinner-validation-comparison.cxx \
                                               ../../src/read_particles_from_TAUOLA.cxx

tauspinner_validation_plots_exe_SOURCES = tauspinner-validation-plots.cxx \
                                          ../../src/read_particles_from_TAUOLA.cxx

test_bornAFB_bornAFB_test_exe_SOURCES = test-bornAFB/bornAFB-test.cxx \
                                        ../../src/read_particles_from_TAUOLA.cxx

test_ipol_test_ipol_exe_SOURCES = test-ipol/test-ipol.cxx \
                                  ../../src/read_particles_from_TAUOLA.cxx

INCLUDES += -I$(PYTHIA8_DIR)/include -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD    += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)


endif

endif

