includedir=$(prefix)/include/TauSpinner

lib_LTLIBRARIES = libTauolaTauSpinner.la

libTauolaTauSpinner_la_SOURCES =            \
src/tau_reweight_lib.cxx                    \
src/nonSM.cxx                               \
src/disth.f                                 \
src/VBF/VBF_init.f      \
src/VBF/VBF_distr.f     \
src/VBF/VBF_functions.f \
src/VBF/vbfdistr.cxx    \
src/VBF/DD.f  \
src/VBF/UU.f  \
src/VBF/UD.f  \
src/VBF/GG.f  \
src/VBF/GU.f  \
src/VBF/GD.f  \
src/VBF/DDX.f \
src/VBF/UUX.f \
src/VBF/UDX.f \
src/VBF/CC.f  \
src/VBF/CD.f  \
src/VBF/CSX.f  \
src/VBF/CUX.f  \
src/VBF/CDX.f  \
src/VBF/DSX.f  \
src/VBF/SD.f  \
src/VBF/SS.f  \
src/VBF/SU.f  \
src/VBF/UCX.f  \
src/VBF/USX.f  \
src/VBF/CCX.f  \
src/VBF/CS.f  \
src/VBF/CU.f  \
src/VBF/DC.f  \
src/VBF/DS.f  \
src/VBF/US.f  \
src/VBF/SCX.f  \
src/VBF/SDX.f  \
src/VBF/SSX.f  \
src/VBF/SUX.f  \
src/VBF/DCX.f

libTauolaTauSpinner_la_FFLAGS = -fno-automatic -fno-backslash -ffixed-line-length-132
libTauolaTauSpinner_la_LIBADD = ../tauola-fortran/libTauolaFortran.la ../src/libTauolaCxxInterface.la

include_HEADERS =                     \
include/TauSpinner/Particle.h         \
include/TauSpinner/SimpleParticle.h   \
include/TauSpinner/Tauola_wrapper.h   \
include/TauSpinner/tau_reweight_lib.h \
include/TauSpinner/nonSM.h            \
include/TauSpinner/vbfdistr.h

INCLUDES = -I$(top_srcdir)/include  \
           -I$(LHAPDF_DIR)/include  \
           -Iinclude

