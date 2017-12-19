HEPMCLOCATION=    /home/erichter/teacher/Dydaktyka2014/TCS-2014/references/TauCP/TauSpinner/hepmc269
MCTESTERLOCATION= /home/erichter/teacher/Dydaktyka2014/TCS-2014/references/TauCP/TauSpinner/MC-TESTER

ifneq ($(MCTESTERLOCATION), )

  ROOTLIB       :=  $(shell root-config --glibs)
  ROOTINC       :=  $(shell root-config --cflags)
  ROOTLIBPATH   :=  $(shell root-config --libdir)
  MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester

  override CFLAGS += -D_USE_MCTESTER_
endif

