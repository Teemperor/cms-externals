//
// Example code which estimates the round-off noise of the FFT code
//
// I. Volobouev
// May 2009

#include <cmath>

// FFTJet headers we will need
#include "fftjet/Grid2d.hh"

// Headers from the "examples" directory
#include "fftjet_typedefs.hh"
#include "SimpleGenerator.hh"
#include "SimpleEvent.hh"
#include "discretizeMCEvent.hh"

int main(int, char**)
{
    // Configure event topology and other job parameters

    // The number of events to process
    const unsigned nevents = 1000;

    // The number of jets to generate per event
    const unsigned njets = 4;

    // Parameters of the jets to generate: minimum and maximum
    // values of pt, eta, phi
    const double minPt = 1.0;
    const double maxPt = 150.0;
    const double minEta = -1.0;
    const double maxEta = 1.0;
    const double minPhi = 0.0;
    const double maxPhi = 2.0*M_PI;

    // Which parton will be used to generate the jets?
    const unsigned pythiaPartonCode = 2;

    // Define the parameters of the "calorimeter" simulation.
    // Note that the eta extent of the calorimeter grid should be
    // wider than the eta extent of the jet spectrum. It should be
    // wide enough so that the "energy leakage" from one side
    // (due to convolutions) does not affect peak finding on the
    // other side.
    const unsigned nEtaCells = 64;
    const unsigned nPhiCells = 64;
    const double caloEtaMax = 3;
    const double caloEtaMin = -caloEtaMax;

    // Where do we place the calorimeter, and what is the magnetic field?
    // The settings below are similar to the conditions found at CMS.
    const double caloRadius = 2.0;
    const double magneticFieldB = 3.8;

    // Simulate calorimeter noise. It will be gaussian, with
    // the standard deviation of the energy given below (in GeV).
    const double caloNoise = 0.15;

    // Calorimeter energy threshold. We will not use cells
    // with energies below this threshold.
    const double caloThreshold = 2.5*caloNoise;

    /////////////////////////////////////////////////////////////////////////

    // Create the MC jet generator
    SimpleGenerator generator(njets, pythiaPartonCode, 
                              minPt, maxPt, minEta, maxEta, minPhi, maxPhi);

    // The structure which will hold the MC event info
    SimpleEvent mcEvent;

    // The structure which will hold the "calorimeter" information
    fftjet::Grid2d<Real> calo(nEtaCells, caloEtaMin, caloEtaMax, nPhiCells, 0.0);

    // Another such structure which will hold doubly transformed image
    fftjet::Grid2d<Real> calo2(calo);

    // Accumulator for the noise data
    fftjet::StatAccumulator noiseAccumulator;

    // The FFT engine
    MyFFTEngine engine(nEtaCells, nPhiCells);

    // Allocate the array which will hold the image data
    Complex* image = engine.allocateComplex();
    assert(image);

    // Cycle over events
    for (unsigned iev=0; iev<nevents; ++iev)
    {
        // Fill the MC event info
        generator.generate(&mcEvent);

        // Discretize the particle energies
        discretizeMCEvent(mcEvent, &calo, magneticFieldB,
                          caloRadius, caloNoise, caloThreshold);

        // Transform the data forward
        engine.transformForward(calo.data(), image);

        // Transform the data back
        engine.transformBack(image, const_cast<Real*>(calo2.data()));

        // Subtract the doubly transformed image from the original
        calo -= calo2;

        // Accumulate the noise statistics
        calo.accumulateDataStats(&noiseAccumulator);
    }

    // Free the image memory
    engine.destroyComplex(image);

    // Print the noise statistics
    const double min = noiseAccumulator.min();
    const double max = noiseAccumulator.max();
    std::cout << "\nFFT noise statistics for " << nevents << " events:\n";
    std::cout << "Min value is " << min << std::endl;
    std::cout << "Max value is " << max << std::endl;
    std::cout << "Mean is " << noiseAccumulator.mean() << std::endl;
    std::cout << "RMS is " << noiseAccumulator.stdev() << std::endl;

    std::cout << "\nPeak finder pulse height cutoff should be at least "
              << 10.0*std::max(fabs(min), fabs(max)) << std::endl;

    return 0;
}
