# This makefile is used to copy files during the "make install" step

if USE_GMAKE
  MAKE = gmake
endif

# files to distribute
EXTRA_DIST = \
    DECAY.DEC \
    isadecay.dat \
    isaparticles.dat \
    evt.pdl

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(prefix)/examples/@PACKAGE@/data
	for file in $(EXTRA_DIST); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(prefix)/examples/@PACKAGE@/data/$$file; done
