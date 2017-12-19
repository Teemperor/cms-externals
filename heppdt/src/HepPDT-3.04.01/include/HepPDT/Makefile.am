## Process this file with automake to produce Makefile.in

if USE_GMAKE
  MAKE = gmake
endif

pkgincludedir = $(prefix)/include/HepPDT
COPY_P = @COPY_P@

pkginclude_HEADERS = 		\
	Constituent.hh		\
	DefTable.hh		\
	HeavyIonUnknownID.hh	\
	Measurement.hh		\
	Measurement.icc		\
	ParticleID.hh		\
	ParticleData.hh		\
	ParticleData.icc	\
	ParticleDataTable.hh	\
	ParticleDataTable.icc	\
	ParticleDataTableComparison.hh	\
	ProcessUnknownID.hh	\
	ResonanceStructure.hh	\
	SimpleProcessUnknownID.hh	\
	SpinState.hh		\
	SpinState.icc		\
	stringtodouble.hh	\
	TableBuilder.hh		\
	TableBuilder.icc	\
	TempParticleData.hh	\
	Version.hh		\
	defs.h


#	CommonParticleData.hh	\
#	CommonParticleData.icc	\
#	DecayChannel.hh		\
#	DecayChannel.icc	\
#	DecayData.hh		\
#	DecayData.icc		\
#	QQChannel.hh		\
#	QQDecay.hh		\
#

# Identify generated file(s) to be removed when 'make clean' is requested:
CLEANFILES = defs.h
