SUBDIRS = generator kinematics

noinst_LTLIBRARIES = libexample.la

libexample_la_SOURCES = discretizeMCEvent.cc processResults.cc

noinst_HEADERS = discretizeMCEvent.hh fftjet_typedefs.hh \
                 FindNClusters.hh iterateRecombinationAlg.hh \
                 matchOneToOne.hh processResults.hh VBuilders.hh \
                 iterateRecombinationAlg.icc

noinst_PROGRAMS = singleScaleExample multiScaleExample \
                  estimateFFTNoise clusteringTreeVis

singleScaleExample_SOURCES = singleScaleExample.cc

multiScaleExample_SOURCES = multiScaleExample.cc FindNClusters.cc

estimateFFTNoise_SOURCES = estimateFFTNoise.cc

clusteringTreeVis_SOURCES = clusteringTreeVis.cc

LDADD = libexample.la ../fftjet/libfftjet.la generator/libgenerator.la \
        kinematics/librk.la $(DEPS_LIBS) $(FLIBS)

AM_CPPFLAGS = -I@top_srcdir@/ $(DEPS_CFLAGS) -Igenerator -Ikinematics

EXTRA_DIST = patreco tables

TESTS = $(noinst_PROGRAMS)
