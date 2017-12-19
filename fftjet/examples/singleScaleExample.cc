//
// Example code which illustrates the usage of the FFTJet package
// in the single-scale jet reconstruction mode. In this simple
// mode FFTJet has two important advantages over, for example,
// the cone algorithm:
//
//  -- The split-merge stage is avoided. Instead, much more
//     useful pattern recognition stage is run before energy
//     recombination is performed.
//
//  -- In the presence of strong magnetic field, jet shapes
//     are notably wider in phi than in eta. FFTJet can take
//     this into account.
//
// I. Volobouev
// April 2009

// FFTJet headers we will need
#include "fftjet/Grid2d.hh"
#include "fftjet/ConstScaleReconstruction.hh"
#include "fftjet/Kernels.hh"
#include "fftjet/PeakSelectors.hh"
#include "fftjet/GaussianNoiseMembershipFcn.hh"
#include "fftjet/KernelRecombinationAlg.hh"
#include "fftjet/KernelConvolver.hh"

// Headers from the "examples" directory
#include "fftjet_typedefs.hh"
#include "SimpleGenerator.hh"
#include "SimpleEvent.hh"
#include "discretizeMCEvent.hh"
#include "processResults.hh"

int main(int, char**)
{
    // Configure event topology and other job parameters

    // The number of events to process
    const unsigned nevents = 20;
    
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
    const double caloEtaMax = M_PI;
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
    //                                                                     //
    //    The parameters below define the FFTJet algorithm behavior        //
    //                                                                     //
    /////////////////////////////////////////////////////////////////////////

    // The ratio between eta and phi bandwidth values for the pattern
    // recognition filter. Naturally, optimal value of this parameter
    // depends on the magnetic field. The value used below is close to
    // optimal for the magnetic field used in this example.
    const double etaToPhiBandwidthRatio = 1.0/3.0;

    // The ratio between bandwidth values used in recombination
    // and in pattern recognition. The optimal ratio will depend
    // on the types of kernel used.
    const double recoToPatBandwidthRatio = 3.0;

    // The scale parameter to use. Bigger scale corresponds to wider
    // kernels used both at the pattern recognition and at the energy
    // recombination stages.
    const double scale = 0.15;

    // Cutoff magnitude for the peaks found at the pattern recognition
    // stage (they should be high enough in order to pass through to
    // recombination). This cutoff should really be optimized depending
    // on the user preference for high efficiency or low noise.
    const double minPeakMagnitude = 0.15;

    /////////////////////////////////////////////////////////////////////////

    // In the single-scale reconstruction mode, the goal of the FFTJet
    // package configuration is to make a ConstScaleReconstruction
    // object which will then take care of running the algorithm.
    // In order to do this, we will need to build a number of lower level
    // objects which specify various aspects of the algorithm behavior.
    // Here is the object hierarchy we are constructing (not including
    // several parameters expressed using built-in types):
    //
    // ConstScaleReconstruction
    //     AbsKernelConvolver (here, KernelConvolver)
    //         AbsFFTEngine (here, MyFFTEngine from "fftjet_typedefs.hh")
    //         AbsKernel2d  (here, Gauss2d)
    //     AbsPeakSelector (here, SimplePeakSelector)
    //     PeakFinder
    //     AbsRecombinationAlg (here, KernelRecombinationAlg)
    //         ScaleSpaceKernel (here, Linear2d)
    //         AbsBgFunctor (from "fftjet_typedefs.hh", here,
    //                       GaussianNoiseMembershipFcn)
    //
    // The classes whose names start with "Abs" are "abstract classes"
    // whose main purpose is to describe functional interfaces required
    // to operate higher level objects. The concrete classes, which
    // implement these interfaces in this particular example, are given in
    // parentheses. Naturally, the implementation classes can be swapped
    // with minimal programming effort. The FFTJet package provides
    // especially numerous implementations of "AbsKernel2d" and
    // "AbsRecombinationAlg" interfaces (together with factories which
    // can be used to create objects of different types from the type name
    // represented by std::string -- consult the user manual for details).
    // User-defined implementation classes can be seamlessly added as well.
    // In other words, FFTJet class system is highly extensible.

    // Now, we will start constructing the relevant objects. First,
    // we will build the pattern recognition kernel (AbsKernel2d instance).
    // It will be a Gaussian with unequal width in eta and phi (the jets
    // are wider in phi due to the presence of magnetic field). The
    // bandwidth values calculated below will be later multiplied by
    // the scale parameter to obtain the actual Gaussian sigmas.
    const double etaBandwidth = sqrt(etaToPhiBandwidthRatio);
    const double phiBandwidth = 1.0/etaBandwidth;

    // Here comes an important fine point. FFT assumes that the grid extent
    // in eta is also 2*Pi, just like in phi. Because of this, when
    // convolutions with the kernel are performed, the eta-phi space
    // will be effectively contracted or expanded in the eta direction.
    // We need to multiply the eta bandwidth by the expansion factor
    // in order to get the requested bandwidth in the original space.
    const double expansionFactor = 2.0*M_PI/(caloEtaMax - caloEtaMin);
    fftjet::Gauss2d kernel(etaBandwidth*expansionFactor, phiBandwidth, 1);

    // Build the FFT engine
    MyFFTEngine engine(nEtaCells, nPhiCells);

    // We can now build the kernel convolver for pattern recognition.
    // The convolver will manage kernel images which do not change
    // from event to event, and will convolve them with data images.
    fftjet::KernelConvolver<Real,Complex> convolver(&engine, &kernel);

    // Build the peak finder. It is important to raise the peak height
    // cutoff above the level of the FFT round-off noise.
    fftjet::PeakFinder peakFinder(1.e-10);

    // Build the peak selector. This selector will filter peaks
    // found by the peak finder before they get to the energy
    // recombination stage. Here, we will use a selector which
    // removes all peaks below a certain magnitude.
    fftjet::SimplePeakSelector selector(minPeakMagnitude);

    // Build the jet membership function for the recombination stage.
    // Here, we will use a very simple kernel which looks, essentially,
    // like an elliptical cone. The constructor looks like this:
    // Linear2d(a, b, p). "a" and "b" are factors for the scale and "p"
    // is the scale power. a*pow(scale,p) gives us the cone semi-major
    // axis in eta, b*pow(scale,p) gives us the semi-major axis in phi.
    fftjet::Linear2d jetMemberFcn(etaBandwidth*recoToPatBandwidthRatio,
                                  phiBandwidth*recoToPatBandwidthRatio, 1);

    // Function which describes "noise cluster" membership.
    // The precise parameters of this function are not important
    // here because the cone membership functions of the jets will
    // almost surely be larger inside the cone and smaller outside.
    fftjet::GaussianNoiseMembershipFcn noiseMemberFcn(1.e-8, 0.0);

    // Build the recombination algorithm. A variety of algorithms
    // can be used. See the user manual for details. We are not
    // going to use fuzzy culstering here: this kind of clustering
    // needs a reasonable probabilistic model for a cell to belong
    // to every jet, and modeling jets by cones is way too crude.
    fftjet::KernelRecombinationAlg<Real,VectorLike,BgData,VBuilder> recoAlg(
        &jetMemberFcn, &noiseMemberFcn, 0.0, 0.0, true, false, false);

    // Build the driver object which runs the FFTJet algorithm sequence
    // for a single scale. Once we have this object, FFTJet is ready to roll.
    fftjet::ConstScaleReconstruction<Real,VectorLike,BgData> sequencer(
        &convolver, &selector, peakFinder, &recoAlg);

    /////////////////////////////////////////////////////////////////////////

    // Create the MC jet generator
    SimpleGenerator generator(njets, pythiaPartonCode, 
                              minPt, maxPt, minEta, maxEta, minPhi, maxPhi);

    // The structure which will hold the MC event info
    SimpleEvent mcEvent;

    // The structure which will hold the "calorimeter" information
    fftjet::Grid2d<Real> calo(nEtaCells, caloEtaMin, caloEtaMax, nPhiCells, 0.0);

    // The vector of reconstructed jets (we will refill it in every event)
    std::vector<fftjet::RecombinedJet<VectorLike> > recoJets;

    // Unclustered energy (vector and scalar) which we will refill
    // in every event. The "unclusScalar" will be the scalar sum of
    // the grid values not assigned to any jet. Due to the way
    // energy discretization is handled in this example, it will be
    // the scalar sum of transverse energies.
    VectorLike unclustered;
    double unclusScalar;

    // Cycle over events
    for (unsigned iev=0; iev<nevents; ++iev)
    {
        // Fill the MC event info
        generator.generate(&mcEvent);

        // Discretize the particle energies
        discretizeMCEvent(mcEvent, &calo, magneticFieldB,
                          caloRadius, caloNoise, caloThreshold);

        // Run the single-scale version of the FFTJet algorithm.
        // If we were to use a more detailed model for jet membership,
        // we would also need to provide a reasonable description
        // of background/noise. In this example, however, the noise
        // information is simply ignored. The cone recombination is
        // sufficiently robust and operates reasonably well in this
        // mode (the price for this robustness is reduction in precision).
        BgData ignored = 0.0;
        int status = sequencer.run(scale, calo, &ignored, 1, 1,
                                   &recoJets, &unclustered, &unclusScalar);
        if (status)
        {
            std::cerr << "Event " << iev
                      << " : ERROR in ConstScaleReconstruction, status "
                      <<  status << std::endl;
            break;
        }

        // In principle, we are done here. We now have the vector
        // of reconstructed jets and the unclustered energy.
        // Run some post-processing steps and print the results.
        processResults(iev, etaToPhiBandwidthRatio, mcEvent,
                       recoJets, unclustered, unclusScalar);
    }

    return 0;
}
