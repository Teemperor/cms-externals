// This code simulated a response of an ideal calorimeter (with some
// readout noise included) in the presence of a magnetic field. Particle
// energies are discretized using pseudorapidity and Et.
//
// I. Volobouev
// May 2009

#ifndef DISCRETIZEMCEVENT_HH_
#define DISCRETIZEMCEVENT_HH_

#include "fftjet/Grid2d.hh"
#include "fftjet_typedefs.hh"

struct SimpleEvent;

// The following function will be used to discretize MC events
void discretizeMCEvent(const SimpleEvent& mcEvent, fftjet::Grid2d<Real>* calo,
                       double magneticFieldB, double calorimeterRadius,
                       double calorimeterNoise, double calorimeterThreshold);

#endif // DISCRETIZEMCEVENT_HH_
