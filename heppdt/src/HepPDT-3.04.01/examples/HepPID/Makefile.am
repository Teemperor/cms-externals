# This makefile is used to copy files during the "make install" step

if USE_GMAKE
  MAKE = gmake
endif

pkgexamdir = $(prefix)/examples/HepPID

# files to distribute
EXTRA_DIST = \
    examListHerwig.cc \
    examListHerwigInterface.F \
    examListHerwig.out.save \
    examListIsajet.cc \
    examListIsajetInterface.F \
    examListIsajet.out.save \
    examListIsajet.lpt.save \
    examListPythia.cc \
    examListPythiaInterface.F \
    examListPythia.out.save \
    examListPythia.lpt.save

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(pkgexamdir)
	$(INSTALL_DATA) GNUmakefile.example $(DESTDIR)$(pkgexamdir)/GNUmakefile
	for file in $(EXTRA_DIST); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(pkgexamdir)/$$file; done
