AC_PREREQ([2.59])

AC_INIT([FFTJet], [1.5.1])

AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIR([m4])

dnl make sure FFTW is installed
PKG_CHECK_MODULES(DEPS, fftw3 >= 3.1.2)
AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_LIBS)

AC_PROG_CXX

dnl options needed to compile Pythia from sources and then link it
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS

AC_PROG_INSTALL

dnl set the default destination directory
AC_PREFIX_DEFAULT([/usr/local])

dnl enable libtool
AC_PROG_LIBTOOL

AC_OUTPUT([Makefile fftjet/Makefile examples/Makefile examples/kinematics/Makefile examples/generator/Makefile fftjet.pc])
