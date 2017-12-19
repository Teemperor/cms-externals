## Process this file with automake to produce Makefile.in

if USE_GMAKE
  MAKE = gmake
endif

INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include

libHepPID_a_SOURCES = \
	ParticleName.cc		\
	translateEvtGen.cc	\
	translateGeanttoPDT.cc	\
	translateHerwig.cc	\
	translateIsajet.cc	\
	translatePythia.cc	\
	translatePDG.cc		\
	translatePDTtoGeant.cc	\
	translateQQ.cc		\
	ParticleIDMethods.cc	\
	Version.cc

libHepPID_so_OBJECTS = $(patsubst %.cc,$(shareddir)/%.$(OBJEXT),$(libHepPID_a_SOURCES))

# these if/else blocks will be parsed by automake and translated
if BUILD_SHARED
  shareddir = shared
if BUILD_VISUAL
  lib_shared = HepPID.$(SHEXT)
else
  lib_shared = libHepPID.$(SHEXT)
endif
  SHFLAGS = @MY_SHFLAGS@
  SHLINK = @MY_SHLINK@
  SHNAME = @MY_SHNAME@
endif
if BUILD_STATIC
  lib_LIBRARIES = libHepPID.a
endif

CLEANFILES =  $(lib_shared)
CYGPATH_W = @CYGPATH_W@

all-local: $(shareddir) $(lib_shared)

$(shareddir):
	test -d $(shareddir) || mkdir $(shareddir)

if BUILD_VISUAL
$(lib_shared): $(libHepPID_so_OBJECTS)
	lib $(SHLINK) $(SHNAME)`$(CYGPATH_W) '$@'` $(libHepPID_so_OBJECTS)

$(shareddir)/%.$(OBJEXT): %.cc
	$(CXXCOMPILE) $(SHFLAGS) /c /Fo$@ `$(CYGPATH_W) '$<'`
else
$(lib_shared): $(libHepPID_so_OBJECTS)
	$(CXXLINK) $(SHLINK) $(SHNAME)$@ $(libHepPID_so_OBJECTS) -o $@

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
