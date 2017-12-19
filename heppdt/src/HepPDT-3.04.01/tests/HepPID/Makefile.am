# Applicable to each test program:

if USE_GMAKE
  MAKE = gmake
endif
INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include

# libraries to link, and whence
if BUILD_STATIC
  LDADD = $(top_builddir)/src/HepPID/libHepPID.a
else
if BUILD_VISUAL
  LDADD = $(top_builddir)/src/HepPID/HepPID.$(SHEXT)
else
  LIBS += $(MY_LD)$(top_builddir)/src/HepPID
  LDADD = $(top_builddir)/src/HepPID/libHepPID.$(SHEXT)
endif
endif

# Identify executables needed during testing:
check_PROGRAMS = \
        listParticleNames testParticleIDMethods \
	listHerwigTranslation listPythiaTranslation listQQTranslation \
	listIsajetTranslation listEvtGenTranslation listPDGTranslation

check_SCRIPTS = \
        listParticleNames.sh testParticleIDMethods.sh  \
	listHerwigTranslation.sh listPythiaTranslation.sh listQQTranslation.sh \
	listIsajetTranslation.sh listEvtGenTranslation.sh listPDGTranslation.sh

# Identify test(s) to run when 'make check' is requested:
TESTS = \
        listParticleNames.sh testParticleIDMethods.sh \
	listHerwigTranslation.sh listPythiaTranslation.sh listQQTranslation.sh \
	listIsajetTranslation.sh listEvtGenTranslation.sh listPDGTranslation.sh

# Identify the test(s) for which failure is the intended outcome:
XFAIL_TESTS = 

# Identify the dependencies on a per-test basis:

listParticleNames_SOURCES = listParticleNames.cc

testParticleIDMethods_SOURCES = testParticleIDMethods.cc

listEvtGenTranslation_SOURCES = listEvtGenTranslation.cc

listHerwigTranslation_SOURCES = listHerwigTranslation.cc

listPythiaTranslation_SOURCES = listPythiaTranslation.cc

listPDGTranslation_SOURCES = listPDGTranslation.cc

listIsajetTranslation_SOURCES = listIsajetTranslation.cc

listQQTranslation_SOURCES = listQQTranslation.cc

# Identify input data file(s) and prototype output file(s):
EXTRA_DIST = \
        listParticleNames.output testParticleIDMethods.output \
	listPythiaTranslation.output listQQTranslation.output \
	listIsajetTranslation.output listEvtGenTranslation.output \
	listHerwigTranslation.output listPDGTranslation.output

# Identify generated file(s) to be removed when 'make clean' is requested:
CLEANFILES = *.out
