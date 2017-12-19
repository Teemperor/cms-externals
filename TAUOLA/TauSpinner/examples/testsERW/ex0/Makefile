include make.inc
include ../../../make.inc

# INFO: Options for program, see README
# USE_SPIN_VALIDATION_ONLY = yes

# INFO: Turn on useful debug printouts
#override CFLAGS+= -DDEBUG_MODE

# INFO: Initialization of option  ' USE_SPIN_VALIDATION_ONLY'
ifdef USE_SPIN_VALIDATION_ONLY
  override CFLAGS += -D_SPIN_VALIDATION_ONLY_
endif

HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
TAUOLA_LIBS    = -L$(PREFIX)/lib -L../lib -lTauolaTauSpinner -lTauolaFortran -lTauolaCxxInterface
LHAPDF_LIBS    = -L$(LHAPDFLOCATION)/lib -lLHAPDF

TAU_SPINNER_HEPMC_INTERFACE = eventReader.o 

all:
ifneq ($(HEPMCLOCATION), )
	@make test_ex0.exe
endif
	@echo "##################################################################"
	@echo " Examples compiled successfully. Available examples:"
	@if test -e test_ex0.exe; then \
	echo " ./test_ex0.exe"; fi
	@echo "(Available examples depend on paths provided during configuration)"
	@echo "##################################################################"

%.o: %.f
	$(F77) $(FCFLAGS) -c $< -o $@

%.o: %.cxx
	$(CC) $(CFLAGS) -c $< -o $@ -I../include $(ROOTINC) \
	-I$(HEPMCLOCATION)/include    \
	-I$(PREFIX)/include           \
	-I$(MCTESTERLOCATION)/include \
	-I$(LHAPDFLOCATION)/include

%.exe: %.o $(TAU_SPINNER_HEPMC_INTERFACE)
	$(CC) -o $@ $^ $(HEPMC_LIBS) $(TAUOLA_LIBS) $(ROOTLIB) $(MC_TESTER_LIBS) $(LHAPDF_LIBS) $(TAUSPINNER_LIBS) \
	-Wl,-rpath,$(PWD)/lib              \
	-Wl,-rpath,$(HEPMCLOCATION)/lib    \
	-Wl,-rpath,$(PREFIX)/lib           \
	-Wl,-rpath,$(ROOTLIBPATH)          \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib \
	-Wl,-rpath,$(LHAPDFLOCATION)/lib \
	-lgfortran

clean:
	rm -f *.o *.exe $(TAU_SPINNER_HEPMC_INTERFACE)

Clean: clean
	rm -f make.inc

make.inc:
	@echo ""
	@echo "Execute ./configure first"
	@echo ""
	@false

../../make.inc:
	@echo ""
	@echo "Compile TauSpinner library first"
	@echo ""
	@false
