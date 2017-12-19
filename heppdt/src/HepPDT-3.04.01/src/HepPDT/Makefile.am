## Process this file with automake to produce Makefile.in

if USE_GMAKE
  MAKE = gmake
endif

INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include

libHepPDT_a_SOURCES = \
	Constituent.cc		\
	DefTable.cc		\
	HeavyIonUnknownID.cc	\
	ParticleID.cc		\
	ParticleDataTable.cc	\
	ProcessUnknownID.cc	\
	ResonanceStructure.cc	\
	TempParticleData.cc	\
	Version.cc		\
	addParticleTable.cc	\
	addPDGParticles.cc	\
	addEvtGenParticles.cc	\
	addIsajetParticles.cc	\
	addHerwigParticles.cc	\
	addPythiaParticles.cc	\
	addQQParticles.cc	\
	calculateWidthFromLifetime.cc	\
	convertTemporaryMap.cc	\
	getIsajetID.cc		\
	getPDGpid.cc		\
	getPythiaid.cc		\
	hasMethods.cc		\
	lifetime.cc		\
	quarks.cc		\
	spindtoi.cc		\
	spinitod.cc		\
	stringtodouble.cc	\
	write.cc

libHepPDT_so_OBJECTS = $(patsubst %.cc,$(shareddir)/%.$(OBJEXT),$(libHepPDT_a_SOURCES))

# these if/else blocks will be parsed by automake and translated
if BUILD_SHARED
  shareddir = shared
if BUILD_VISUAL
  lib_shared = HepPDT.$(SHEXT)
else
  lib_shared = libHepPDT.$(SHEXT)
  LDADD = -L${top_builddir}/src/HepPID -lHepPID
endif
  SHFLAGS = @MY_SHFLAGS@
  SHLINK = @MY_SHLINK@
  SHNAME = @MY_SHNAME@
endif
if BUILD_STATIC
  lib_LIBRARIES = libHepPDT.a
endif

#libHepPDT_a_LIBADD = -L${top_builddir}/src/HepPID -lHepPID

CLEANFILES =  $(lib_shared)
CYGPATH_W = @CYGPATH_W@

all-local: $(shareddir) $(lib_shared)

$(shareddir):
	test -d $(shareddir) || mkdir $(shareddir)

if BUILD_VISUAL
$(lib_shared): $(libHepPDT_so_OBJECTS)
	lib $(SHLINK) $(SHNAME)`$(CYGPATH_W) '$@'` $(libHepPDT_so_OBJECTS)

$(shareddir)/%.$(OBJEXT): %.cc
	$(CXXCOMPILE) $(SHFLAGS) /c /Fo$@ `$(CYGPATH_W) '$<'`
else
$(lib_shared): $(libHepPDT_so_OBJECTS)
	$(CXXLINK) $(SHLINK) $(SHNAME)$@ $(libHepPDT_so_OBJECTS) $(LDADD) -o $@

$(shareddir)/%.$(OBJEXT): %.cc
	$(CXXCOMPILE) $(SHFLAGS) -c -o $@ $<
endif

install-exec-local: $(lib_shared)
	@$(NORMAL_INSTALL)
	test -z "$(libdir)" || $(mkdir_p) "$(DESTDIR)$(libdir)"
	@list='$(lib_shared)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo " $(INSTALL) '$$p' '$(DESTDIR)$(libdir)/$$p'"; \
	    $(INSTALL) "$$p" "$(DESTDIR)$(libdir)/$$f"; \
	  else :; fi; \
	done

clean-local:
	-test -d $(shareddir) && rm -f $(shareddir)/*.$(OBJEXT)

distclean-local:
	-test -d $(shareddir) && rm -rf ./$(shareddir)
