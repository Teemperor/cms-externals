# ======================================================================
# This is the HepPDT configuration input file
# Process this file with autoconf to produce a configure script.
# ======================================================================


# ----------------------------------------------------------------------
# Identify the package and initialize the autotools:
# ----------------------------------------------------------------------

AC_PREREQ(2.59)
AC_INIT([HepPDT],[3.04.01],[https://savannah.cern.ch/projects/heppdt/],[HepPDT])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE(1.9 foreign)
AC_CONFIG_SRCDIR([src/HepPDT/Constituent.cc])
AC_CONFIG_SRCDIR([src/HepPID/ParticleName.cc])
AC_CONFIG_HEADER([include/HepPDT/defs.h])

# ----------------------------------------------------------------------
# --enable-gmake
# ----------------------------------------------------------------------
AC_ARG_ENABLE(gmake,
   AC_HELP_STRING([--enable-gmake],[use gmake (default is make)]),
        [case "${enableval}" in
          yes) use_gmake="yes";;
          no)  use_gmake="no";;
          *)   use_gmake="no";;
        esac],
        [use_gmake="no"])

# ----------------------------------------------------------------------
# define --enable-visual
# ----------------------------------------------------------------------
AC_ARG_ENABLE(visual,
   AC_HELP_STRING([--enable-visual],[on by default when using Visual C++]),
        [case "${enableval}" in
          yes) build_visual="yes";;
          no)  build_visual="no";;
          *)   build_visual="no";;
        esac],
        [build_visual="no"])

# ----------------------------------------------------------------------
# define --enable-shared and --enable-static
# ----------------------------------------------------------------------
AC_ARG_ENABLE(shared,
   AC_HELP_STRING([--disable-shared],[do not build shared libraries]),
        [case "${enableval}" in
          yes) build_shared="yes";;
          no)  build_shared="no";;
          *)   build_shared="yes";;
        esac],
        [build_shared="yes"])
AC_ARG_ENABLE(static,
   AC_HELP_STRING([--disable-static],[do not build static libraries]),
        [case "${enableval}" in
          yes) build_static="yes";;
          no)  build_static="no";;
          *)   build_static="yes";;
        esac],
        [build_static="yes"])

# ----------------------------------------------------------------------
# Check for needed programs:
# ----------------------------------------------------------------------
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# ----------------------------------------------------------------------
# Checks for libraries.
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Checks for header files.
# ----------------------------------------------------------------------
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([stdint.h])

# ----------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
# ----------------------------------------------------------------------
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_CHECK_TYPES([ptrdiff_t])

# ----------------------------------------------------------------------
# Checks for library functions.
# ----------------------------------------------------------------------
AC_CHECK_FUNCS([pow sqrt])

# ----------------------------------------------------------------------
# OS dependent differences:
# ----------------------------------------------------------------------

# copy and diff
case "$target" in
*-*-win32*)
  COPY_P="copy -p";DIFF_Q="diff -q -b"
  ;;
*-*-cygwin*)
  COPY_P="copy -p";DIFF_Q="diff -q -b"
  ;;
*-*-solaris*)
  COPY_P="cp -p";DIFF_Q="diff -b"
  use_gmake="yes"
  ;;
*-*-linux*)
  COPY_P="cp -p";DIFF_Q="diff -q -b"
  ;;
*-*-darwin*) 
  COPY_P="cp -p";DIFF_Q="diff -q -b"
  ;;
*)
  COPY_P="cp -p";DIFF_Q="diff -q -b"
esac

AM_CONDITIONAL(USE_GMAKE, test x$use_gmake = xyes)

AC_SUBST(COPY_P)
AC_SUBST(DIFF_Q)


# ----------------------------------------------------------------------
# compiler dependent differences
# notice that the default is presumed to be some flavor of g++
# ----------------------------------------------------------------------
case "$CXX" in
cl) 
   AM_CXXFLAGS="-EHsc -nologo -GR -MD"
   CXXFLAGS=" "
   AC_SUBST(CXXFLAGS)
   MY_SHFLAGS=" "
   MY_SHLINK="/NOLOGO";
   MY_SHNAME="/OUT:"
   SHEXT=lib
   build_static="no"
   build_shared="yes"
   build_visual="yes"
   ;;
CC)  
   AM_CXXFLAGS="-O"
   MY_SHFLAGS="-KPIC -DPIC"
   MY_SHLINK="-G"
   MY_SHNAME="-h"
   SHEXT=so
   build_static="yes"
   build_shared="no"
   ;;
*)  
   AM_CXXFLAGS="-O -ansi -pedantic -Wall"
   MY_SHFLAGS="-fPIC -DPIC"
   case "$target" in
   *-*-linux*)
     MY_SHLINK="-shared"
     MY_SHNAME="-Wl,-soname,"
     SHEXT=so
     ;;
   *-*-darwin*) 
     MY_SHLINK="-dynamiclib -single_module"
     MY_SHNAME="-Wl,-install_name,"
     SHEXT=dylib
     ;;
   *)
     MY_SHLINK="-shared"
     MY_SHNAME="-Wl,-soname,"
     SHEXT=so
   esac
esac

AC_SUBST(AM_CXXFLAGS)
AC_SUBST(MY_SHFLAGS)
AC_SUBST(MY_SHLINK)
AC_SUBST(MY_SHNAME)
AC_SUBST(SHEXT)

AM_CONDITIONAL(BUILD_SHARED, test x$build_shared = xyes)
AM_CONDITIONAL(BUILD_STATIC, test x$build_static = xyes)
AM_CONDITIONAL(BUILD_VISUAL, test x$build_visual = xyes)

# ----------------------------------------------------------------------
# Identify the files that 'configure' is to produce:
# ----------------------------------------------------------------------

# Makefiles:
AC_CONFIG_FILES([Makefile
                 include/HepPDT/Makefile
                 include/HepPID/Makefile
                 src/HepPDT/Makefile
                 src/HepPID/Makefile
                 doc/HepPDT/Makefile
                 doc/HepPID/Makefile
                 data/Makefile
                 tests/HepPID/Makefile
                 tests/HepPDT/Makefile
		 tests/HepPDT/testHepPDT.input
		 tests/HepPDT/listPythiaNames.cc
		 tests/HepPDT/listPDGNames.cc
		 tests/HepPDT/listEvtGenNames.cc
		 tests/HepPDT/testReadParticleTable.cc
		 tests/HepPDT/testReadEvtGen.cc
		 tests/HepPDT/testReadIsajet.cc
		 tests/HepPDT/testReadQQ.cc
                 examples/data/Makefile
                 examples/HepPID/Makefile
                 examples/HepPID/GNUmakefile.example
                 examples/HepPDT/Makefile
                 examples/HepPDT/GNUmakefile.example])

AC_CONFIG_FILES([tests/HepPDT/testPID.sh], [chmod +x tests/HepPDT/testPID.sh])
AC_CONFIG_FILES([tests/HepPDT/testHepPDT.sh], [chmod +x tests/HepPDT/testHepPDT.sh])
AC_CONFIG_FILES([tests/HepPDT/listPythiaNames.sh], [chmod +x tests/HepPDT/listPythiaNames.sh])
AC_CONFIG_FILES([tests/HepPDT/listPDGNames.sh], [chmod +x tests/HepPDT/listPDGNames.sh])
AC_CONFIG_FILES([tests/HepPDT/listEvtGenNames.sh], [chmod +x tests/HepPDT/listEvtGenNames.sh])
AC_CONFIG_FILES([tests/HepPDT/testReadEvtGen.sh], [chmod +x tests/HepPDT/testReadEvtGen.sh])
AC_CONFIG_FILES([tests/HepPDT/testReadParticleTable.sh], [chmod +x tests/HepPDT/testReadParticleTable.sh])
AC_CONFIG_FILES([tests/HepPDT/testReadIsajet.sh], [chmod +x tests/HepPDT/testReadIsajet.sh])
AC_CONFIG_FILES([tests/HepPDT/testReadQQ.sh], [chmod +x tests/HepPDT/testReadQQ.sh])
AC_CONFIG_FILES([tests/HepPID/testParticleIDMethods.sh], [chmod +x tests/HepPID/testParticleIDMethods.sh])
AC_CONFIG_FILES([tests/HepPID/listParticleNames.sh], [chmod +x tests/HepPID/listParticleNames.sh])
AC_CONFIG_FILES([tests/HepPID/listIsajetTranslation.sh], [chmod +x tests/HepPID/listIsajetTranslation.sh])
AC_CONFIG_FILES([tests/HepPID/listHerwigTranslation.sh], [chmod +x tests/HepPID/listHerwigTranslation.sh])
AC_CONFIG_FILES([tests/HepPID/listEvtGenTranslation.sh], [chmod +x tests/HepPID/listEvtGenTranslation.sh])
AC_CONFIG_FILES([tests/HepPID/listPythiaTranslation.sh], [chmod +x tests/HepPID/listPythiaTranslation.sh])
AC_CONFIG_FILES([tests/HepPID/listPDGTranslation.sh], [chmod +x tests/HepPID/listPDGTranslation.sh])
AC_CONFIG_FILES([tests/HepPID/listQQTranslation.sh], [chmod +x tests/HepPID/listQQTranslation.sh])

# ----------------------------------------------------------------------
# Finish up:
# ----------------------------------------------------------------------

AC_OUTPUT


# ======================================================================
