# Applicable to each test program:

if USE_GMAKE
  MAKE = gmake
endif
INCLUDES = -I$(top_builddir)/include -I$(top_srcdir)/include

# libraries to link, and whence
if BUILD_STATIC
  LDADD = $(top_builddir)/src/HepPDT/libHepPDT.a \
          $(top_builddir)/src/HepPID/libHepPID.a
else
if BUILD_VISUAL
  LDADD = $(top_builddir)/src/HepPDT/HepPDT.$(SHEXT) \
          $(top_builddir)/src/HepPID/HepPID.$(SHEXT)
else
  LIBS += $(MY_LD)$(top_builddir)/src/HepPDT \
          $(MY_LD)$(top_builddir)/src/HepPID
  LDADD = $(top_builddir)/src/HepPDT/libHepPDT.$(SHEXT) \
          $(top_builddir)/src/HepPID/libHepPID.$(SHEXT)
endif
endif

# Identify executables needed during testing:
check_PROGRAMS = \
        testPID testHepPDT   \
	testReadEvtGen testReadQQ testReadIsajet \
	listPythiaNames listPDGNames listEvtGenNames \
	testReadParticleTable

check_SCRIPTS = \
        testPID.sh testHepPDT.sh testHepPDT.input \
	testReadEvtGen.sh testReadQQ.sh testReadIsajet.sh \
	listPythiaNames.sh listPDGNames.sh listEvtGenNames.sh \
	testReadParticleTable.sh

# Identify test(s) to run when 'make check' is requested:
TESTS = \
        testPID.sh testHepPDT.sh \
	testReadEvtGen.sh testReadQQ.sh testReadIsajet.sh \
	listPythiaNames.sh listPDGNames.sh listEvtGenNames.sh \
	testReadParticleTable.sh

# Identify the test(s) for which failure is the intended outcome:
XFAIL_TESTS = 

# Identify the dependencies on a per-test basis:
testPID_SOURCES = testPID.cc

testHepPDT_SOURCES = testHepPDT.cc

testReadEvtGen_SOURCES = testReadEvtGen.cc

testReadParticleTable_SOURCES = testReadParticleTable.cc

testReadIsajet_SOURCES = testReadIsajet.cc

testReadQQ_SOURCES = testReadQQ.cc

listPythiaNames_SOURCES = listPythiaNames.cc

listPDGNames_SOURCES = listPDGNames.cc

listEvtGenNames_SOURCES = listEvtGenNames.cc

# Identify input data file(s) and prototype output file(s):
EXTRA_DIST = \
        testPID.output testHepPDT.output testHepPDT.input \
	testHepPDTfragment.output testHepPDTstatus.output \
	testHepPDTtable.output TestNuclearFragment.hh \
	testReadEvtGen.output testReadQQ.output testReadIsajet.output \
	testReadParticleTable.output testReadParticleTableStatus.output \
	listQQ.dec listPDGNames.output \
	listPythia.tbl listPythiaNames.output \
	listEvtGenNames.output

# Identify generated file(s) to be removed when 'make clean' is requested:
CLEANFILES = *.out
