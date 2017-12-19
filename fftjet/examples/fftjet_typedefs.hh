// Basic choices which must be made in order to use the FFTJet package

#ifndef FFTJET_TYPEDEFS_HH_
#define FFTJET_TYPEDEFS_HH_

// The header file for the FFTW library
#include "fftw3.h"

// The header file for the relativistic kinematics library
#include "rk.hh"

// Classes which build 4-momenta out of energy and direction
#include "VBuilders.hh"

// Header file for the concrete FFT engine used
#include "fftjet/FFTWDoubleEngine.hh"

// Header file for the functor interface
#include "fftjet/SimpleFunctors.hh"

// The following three typedefs reflect the choice of the
// double precision FFTW library for performing DFFTs
typedef double Real;
typedef fftw_complex Complex;
typedef fftjet::FFTWDoubleEngine MyFFTEngine;

// The next typedef reflects the choice of the 4-vector class
typedef rk::P4 VectorLike;

// The following typedef tells how the 4-vectors will be
// constructed from the grid points in the eta-phi space
typedef PtEtaP4Builder VBuilder;

// The following typedef tells which type (or class) is used
// to provide pileup/background information
typedef double BgData;

// The following typedef defines the interface for the functor
// which calculates the pileup/noise membership function
typedef fftjet::Functor2<double,double,BgData> AbsBgFunctor;

#endif // FFTJET_TYPEDEFS_HH_
