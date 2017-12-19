## Process this file with automake to produce Makefile.in

if USE_GMAKE
  MAKE = gmake
endif

PDG_TABLES = PDG_mass_width_2004.mc \
             mass_width_2006.csv \
	     mass_width_2006.mc \
             mass_width_2008.csv \
	     mass_width_2008.mc \
	     particle.tbl

EXTRA_DIST = $(PDG_TABLES)
      

install-data-local:
	$(mkinstalldirs) $(DESTDIR)$(prefix)/data
	for file in $(PDG_TABLES); do \
	  $(INSTALL_DATA) $(srcdir)/$$file $(DESTDIR)$(prefix)/data/$$file; done

