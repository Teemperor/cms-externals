include ../make.inc
include ../../../make.inc

HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
TAUOLA_LIBS    = -L$(PREFIX)/lib -lTauolaTauSpinner -lTauolaFortran -lTauolaCxxInterface
LHAPDF_LIBS    = -L$(LHAPDFLOCATION)/lib -lLHAPDF

TAU_SPINNER_HEPMC_INTERFACE = ../read_particles_from_TAUOLA.o

all:
ifneq ($(HEPMCLOCATION), )
	@make hepmc-tauola-redecay.exe
	@make applications-plots.exe
	@make applications-comparison.exe
	@make applications-fits/applications-fits.exe
endif

%.exe: %.cxx $(TAU_SPINNER_HEPMC_INTERFACE)
	$(CC) $(CFLAGS) $< -o $@ -I.. -I$(shell root-config --incdir) \
	-I$(HEPMCLOCATION)/include    \
	-I$(PREFIX)/include	   \
	-I$(LHAPDFLOCATION)/include   \
	$(TAU_SPINNER_HEPMC_INTERFACE) $(HEPMC_LIBS) $(TAUOLA_LIBS) $(shell root-config --glibs) $(LHAPDF_LIBS) $(TAUSPINNER_LIBS) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib	  \
	-Wl,-rpath,$(PREFIX)/lib		 \
	-Wl,-rpath,$(shell root-config --libdir) \
	-Wl,-rpath,$(LHAPDFLOCATION)/lib \
	-lgfortran

applications-fits/applications-fits.exe: applications-fits/applications-fits.cxx applications-fits/common.o $(TAU_SPINNER_HEPMC_INTERFACE)
	$(CC) $(CFLAGS) $< applications-fits/common.o -o $@ -I.. -I$(shell root-config --incdir) \
	-I$(HEPMCLOCATION)/include    \
	-I$(PREFIX)/include	   \
	-I$(LHAPDFLOCATION)/include   \
	$(TAU_SPINNER_HEPMC_INTERFACE) $(HEPMC_LIBS) $(TAUOLA_LIBS) $(shell root-config --glibs) $(LHAPDF_LIBS) $(TAUSPINNER_LIBS) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib	  \
	-Wl,-rpath,$(PREFIX)/lib		 \
	-Wl,-rpath,$(shell root-config --libdir) -lRooFit -lRooFitCore -lMinuit \
	-Wl,-rpath,$(LHAPDFLOCATION)/lib \
	-lgfortran

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@ -I.. $(ROOTINC) \
	-I$(HEPMCLOCATION)/include    \
	-I$(PREFIX)/include           \
	-I$(MCTESTERLOCATION)/include \
	-I$(LHAPDFLOCATION)/include

clean:
	rm -f *.o *.exe applications-fits/*.exe
