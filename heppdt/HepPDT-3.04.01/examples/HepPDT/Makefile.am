# This makefile is used to copy files during the "make install" step

if USE_GMAKE
  MAKE = gmake
endif

# files to distribute
EXTRA_DIST = \
    examMyPDT.cc \
    examMyPDT.out.save

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(prefix)/examples/@PACKAGE@
	$(INSTALL_DATA) GNUmakefile.example $(DESTDIR)$(prefix)/examples/@PACKAGE@/GNUmakefile
	for file in $(EXTRA_DIST); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(prefix)/examples/@PACKAGE@/$$file; done
