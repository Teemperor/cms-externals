//
// This example illustrates the usage of the FFTJet package
// in the multiscale jet reconstruction mode. This mode allows
// the users to devise and apply various jet pattern recognition
// strategies in the scale space. Here, a simple pattern recognition
// method is used in which an event is clustered into a certain
// user-defined number of jets. This example also shows how to run
// the energy recombination stage iteratively.
//
// The example uses a rather sophisticated jet membership function
// which is, essentially, the jet fragmentation function with respect
// to eta, phi, and Et variables observed in the calorimeter cells.
// This fragmentation function depends on the energy recombination
// scale parameter taken to be inverse of the jet Pt. The function
// is normalized as follows: its integral over eta, phi, and Et
// equals the average number of calorimeter cells into which a single
// jet deposits its energy in one event. At the same time, the
// integral of this function multiplied by Et equals the average
// deposited jet Et for a particular energy recombination scale.
//
// Unfortunately, the relationship between this "detector-level"
// fragmentation and the theoretical fragmentation functions for
// the jet longitudinal and transverse momenta is far from trivial.
// The "detector-level" fragmentation function includes resolution
// effects, charge particle deviations in the magnetic field, dE/dx
// energy losses, photon conversions, etc. The only general way to
// build such a function in practice is to fully simulate (or observe
// in the experiment) the detector response to a large number of
// well-isolated jets and then average the results. If one knows
// in good detail the shape of the jets one expects to see, one should
// be able to utilize this knowledge in order to improve the jet energy
// reconstruction precision. FFTJet allows you to exploit this knowledge.
//
// Of course, the level of sophistication of jet shape modeling
// can be reduced or increased. For example, eta-phi jet energy
// profiles appear to work well as jet membership functions. They are
// significantly less complex then the function used here (we need
// one less dimension to represent them), and they result in a more
// robust clustering. On the other hand, for ultimate precision,
// one can include more variables into the "detector-level"
// fragmentation function, such as charged particle multiplicity,
// charged energy fraction, assumption about the jet flavor, etc.
// To use such additional information, one has to associate each
// precluster with its own membership function. This association
// is supported by FFTJet design, but the technique is not used
// in this example.
//
// With the pattern recognition strategy employed in this example
// (cluster event into a certain predefined number of jets), it makes
// sense to iterate the energy recombination stage. The jet energies
// and directions at the next iteration are taken from the previous
// iteration. If the jets do not stabilize, it serves as an indicator
// that either the number of jets or the initial configuration were
// chosen poorly. Of course, one has to weight the improvement achieved
// through such iterations against the CPU time needed. The
// computational complexity of every such iteration is O(J*N), where
// J is the number of jets and N is the number of cells in the energy
// discretization grid.
//
// To run the iterations in a meaningful way, we need to know how
// to correct the jets: that is, how to relate the jet energy
// found by the reconstruction algorithm during the previous
// iteration to the actual energy of the jet which determines the
// starting point for the next iteration. When energy losses in the
// material are present or when the magnetic field does not allow
// low-momentum charged particles and gamma conversions to reach the
// calorimeter, this correction becomes essential (and, of course, all
// HEP experiments apply such a correction to all jet reconstruction
// algorithms anyway). This correction can be determined, for example,
// by studying results produced by the algorithm on well-separated jets.
// In this example, the correction will be built on-the-fly from the
// jet response data. "Jet response" is some kind of a statistical
// location measure (mean, median, etc) of the ratio between
// reconstructed and actual jet Pt as a function of actual jet Pt.
//
// In order to declare that the iterations have converged, we will
// define a distance (in the chi-squared sense) between jet momenta
// obtained in two consecutive energy recombination iterations.
// If this distance is small enough, the iterations will be deemed
// converged. The squared distance between jet configurations will be
// defined as the average squared distance between individual jets.
// For two jets, the distance will be defined as
// ((delta eta)/(eta bandwidth))^2 + ((delta phi)/(phi bandwidth))^2 +
// ((delta Pt)/(Pt * (relative Pt bandwidth)))^2. The ratio between 
// eta and phi bandwidth values will be set to "etaToPhiBandwidthRatio",
// and their product will be set to 1. The relative Pt bandwidth
// will be an additional parameter.
//
// I. Volobouev
// May 2009

#include <fstream>

// FFTJet headers we will need
#include "fftjet/Grid2d.hh"
#include "fftjet/ClusteringSequencer.hh"
#include "fftjet/DiscreteGauss2d.hh"
#include "fftjet/PeakSelectors.hh"
#include "fftjet/InterpolatedMembershipFcn.hh"
#include "fftjet/GaussianNoiseMembershipFcn.hh"
#include "fftjet/KernelRecombinationAlg.hh"
#include "fftjet/FrequencyKernelConvolver.hh"
#include "fftjet/EquidistantSequence.hh"
#include "fftjet/ProximityClusteringTree.hh"
#include "fftjet/PeakEtaPhiDistance.hh"
#include "fftjet/LinearInterpolator1d.hh"
#include "fftjet/JetMagnitudeMapper.hh"

// Headers from the "examples" directory
#include "fftjet_typedefs.hh"
#include "SimpleGenerator.hh"
#include "SimpleEvent.hh"
#include "discretizeMCEvent.hh"
#include "processResults.hh"
#include "FindNClusters.hh"
#include "iterateRecombinationAlg.hh"

int main(int, char**)
{
    // Some useful typedefs which simplify subsequent typing
    typedef fftjet::RecombinedJet<VectorLike> RecoJet;
    typedef fftjet::ProximityClusteringTree<fftjet::Peak,long> ClusteringTree;

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

    // How many scales to use, and the range of scales. At each scale,
    // the width of the pattern recognition filter (convolved with the
    // event energy depoisition distribution) will be proportional
    // to the scale.
    const unsigned nScales = 50;
    const double minScale = 0.1;
    const double maxScale = 0.8;

    // Parameters of the scale-dependent noise suppression formula.
    // All peaks found by the peak finder will be suppressed
    // if their magnitude is smaller than a*pow(scale,p) + b.
    // In order to find good values of parameters a, p, and b,
    // one has to examine scale-dependent distributions of peak
    // magnitudes produced by pure noise (0 jets). The parameters
    // should be chosen in such a way that none or just a very small
    // fraction of pure noise peaks pass this selection criterion.
    // Of course, this kind of tuning can be easily performed with
    // any reasonable histogramming and fitting package. Explaining
    // how to do this is beyond the scope of this example.
    const double noise_a = 0.00540657735;
    const double noise_p = -1.49002109;
    const double noise_b = 0.0022395466;

    // Are we going to use crisp or fuzzy clustering at the energy
    // recombination stage? It may be useful to use fuzzy clustering
    // when one expects that a significant fraction of energy deposits
    // is indeed shared between two or more jets (in multi-jet events,
    // when significant smearing is present due to shower development
    // in the calorimeter, etc). In the example given here the use of
    // fuzzy clustering is borderline -- both crisp and fuzzy clustering
    // result in about the same final jet energy determination precision.
    //
    const bool isCrisp = false;

    // The parameter below is the name of the file where the jet
    // membership function is stored. The function is represented
    // by a table of values on a 4-dimensional rectangular grid
    // (stored in a sparse manner). Between the grid points, the
    // function values are interpolated linearly in eta, phi, Et,
    // and log(scale).
    //
    const char* jetMembershipFcnFileName = "tables/pythia_diffshapes_0.04.dat";

    // From time to time we will see the cases when more than one
    // jet contributes energy to the same discretization cell. When
    // 4-dimensional jet membership functions are used, two jets
    // may sometimes create a deposit whose probability will be 0
    // for each of these jets taken alone. FFTJet can detect such
    // cases by finding cells with 0 signal probability and very low
    // noise probability. One can turn on an alternative handling
    // of such cases. In this alternative mode FFTJet will attempt
    // to distributes the energy between two closest jets, up to
    // a certain limit which depends on the jet shape. To turn such
    // handling on, the following parameter should be set to some
    // positive value.
    const double unlikelyNoiseProbability = 1.e-5;

    // Parameters for iterating the jet energy recombination algorithm.
    // To disable the iterations, simply set the "maxIterations" parameter
    // to 0 or 1.
    const char* jetResponseDataFile = "tables/fuzzy_4d_response.dat";
    const unsigned maxIterations = 20;
    const double relativePtBandwidth = 1.0;
    const double convergedDistanceSquared = 2.5e-05;

    // We will also need to provide an initial guess of the jet momentum
    // from the peak magnitude found at the pattern recognition stage.
    // The name of the file which describes the response of the pattern
    // recognition peaks is given here.
    const char* peakResponseDataFile = "tables/peak_response.dat";

    /////////////////////////////////////////////////////////////////////////

    // In the multiscale reconstruction mode, the FFTJet package is
    // configured, essentially, by creating three high level objects
    // with types "ClusteringSequencer", "AbsPatternRecognitionAlg",
    // and "AbsRecombinationAlg". During event processing, the
    // ClusteringSequencer object will be growing "clustering trees",
    // one tree per event. The clustering tree is, basically, a hierarchial
    // clustering dendrogram in which every node is associated with
    // a point in the eta-phi-scale space. The clustering tree conforms
    // to the "AbsClusteringTree" interface, and it is used to perform
    // jet pattern recognition in the scale space. I suggest that you
    // implement your pattern recognition strategy inside a separate
    // class derived from "AbsPatternRecognitionAlg". The preclusters
    // found during the pattern recognition stage are then fed into
    // an energy recombination algorithm derived from "AbsRecombinationAlg"
    // (FFTJet provides a variety of "AbsRecombinationAlg" implementations).
    //
    // To make a ClusteringSequencer, we will need to build several
    // lower level objects which specify various aspects of the algorithm
    // behavior. Here is the object hierarchy we are constructing (not
    // including several parameters expressed using built-in types):
    //
    // ClusteringSequencer
    //     AbsKernelConvolver (here, KernelConvolver)
    //         AbsFFTEngine (here, MyFFTEngine from "fftjet_typedefs.hh")
    //         AbsKernel2d  (here, Gauss2d)
    //     AbsPeakSelector (here, ScalePowerPeakSelector)
    //     PeakFinder
    //     std::vector<double> (here, EquidistantInLogSpace)
    //
    // The object hierarchy for the recombination algorithm looks like this:
    //
    // AbsRecombinationAlg (here, KernelRecombinationAlg)
    //     ScaleSpaceKernel (here, InterpolatedMembershipFcn)
    //     AbsBgFunctor (from "fftjet_typedefs.hh", here,
    //                   GaussianNoiseMembershipFcn)
    //
    // At this time, the FFTJet package does not provide any implementations
    // of "AbsPatternRecognitionAlg" (apart from clustering event energy into
    // a predefined number of jets done in this example). The main reason
    // for this omission is that optimal pattern recognition strategies are
    // going to be experiment and analysis-dependent. General principles for
    // building jet pattern recognition strategies applicable in a wide
    // variety of situations are not obvious, and will be subject of future
    // research (in particular, standard blob detectors in the Gaussian
    // scale-space based on scale-normalized Laplacian and Hessian do not
    // work well for jets because angular jet energy profiles do not decay
    // quickly enough).

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
    fftjet::DiscreteGauss2d kernel(etaBandwidth*expansionFactor, phiBandwidth,
                                   nEtaCells, nPhiCells);

    // Build the FFT engine
    MyFFTEngine engine(nEtaCells, nPhiCells);

    // We can now build the kernel convolver for pattern recognition.
    // The convolver will manage kernel images which do not change
    // from event to event, and will convolve them with data images.
    fftjet::FrequencyKernelConvolver<Real,Complex> convolver(&engine, &kernel);

    // Build the peak finder. It is important to raise the peak height
    // cutoff above the level of the FFT round-off noise.
    fftjet::PeakFinder peakFinder(1.e-10);

    // Build the peak selector. This selector will filter peaks found
    // by the peak finder before they get into the clustering tree
    // (so that they will never reach the energy recombination stage).
    fftjet::ScalePowerPeakSelector selector(noise_a, noise_p, noise_b);

    // Build the set of scales for which the clustering tree
    // will be constructed
    fftjet::EquidistantInLogSpace clusterScales(minScale, maxScale, nScales);

    // Now we can construct the clustering sequencer. In this example,
    // we will cluster events using a predefined set of scales. It may
    // also be interesting to perform adaptive clustering -- see the
    // user manual for details.
    fftjet::ClusteringSequencer<Real> clusteringSequencer(
        &convolver, &selector, peakFinder, clusterScales);

    // Read in the jet membership function for the recombination stage
    fftjet::ScaleSpaceKernel* jetMemFcn = 0;
    {    
        std::ifstream in(jetMembershipFcnFileName, std::ios_base::in |
                                                   std::ios_base::binary);
        if (in.is_open())
            jetMemFcn = fftjet::InterpolatedMembershipFcn<float>::read(in);
    }
    if (jetMemFcn == 0)
    {
        std::cerr << "Failed to read jet membership function from file \""
                  << jetMembershipFcnFileName << '"' << std::endl;
        return 1;
    }

    // Function which describes "noise cluster" membership.
    // If the jet membership function is taken to be the detector-level
    // jet fragmentation function then this should be the noise occupancy
    // per unit energy and per cell so that it has compatible normalization.
    const double cellArea = (caloEtaMax - caloEtaMin)/nEtaCells*
                             2*M_PI/nPhiCells;
    fftjet::GaussianNoiseMembershipFcn noiseMemFcn(1.e-10, 1.0/cellArea);

    // Build the recombination algorithm. A variety of algorithms
    // can be used. See the user manual for details.
    fftjet::KernelRecombinationAlg<Real,VectorLike,BgData,VBuilder> recoAlg(
        jetMemFcn, &noiseMemFcn, unlikelyNoiseProbability,
        0.0, isCrisp, false, false);
    
    // Read in the jet response data and build the correction for
    // the iterative energy recombination
    fftjet::LinearInterpolator1d* jetResponse = 0;
    {
        std::ifstream in(jetResponseDataFile, std::ios_base::in |
                                              std::ios_base::binary);
        if (in.is_open())
            jetResponse = fftjet::LinearInterpolator1d::read(in);
    }
    if (jetResponse == 0)
    {
        std::cerr << "Failed to read jet response data from file \""
                  << jetResponseDataFile << '"' << std::endl;
        return 1;
    }
    fftjet::JetMagnitudeMapper<RecoJet> jetCorrector(
        *jetResponse, jetResponse->xMax(), 1000);

    // Build the pattern recognition algorithm
    const double integFactor = nEtaCells*nPhiCells/(caloEtaMax - caloEtaMin);
    FindNClusters patternRecoAlg(njets, integFactor, peakResponseDataFile);
    if (!patternRecoAlg.isValid())
    {
        std::cerr << "Failed to initialize pattern recognition using file \""
                  << peakResponseDataFile << '"' << std::endl;
        return 1;
    }

    /////////////////////////////////////////////////////////////////////////

    // Create the MC jet generator
    SimpleGenerator generator(njets, pythiaPartonCode, 
                              minPt, maxPt, minEta, maxEta, minPhi, maxPhi);

    // The structure which will hold the MC event info
    SimpleEvent mcEvent;

    // The structure which will hold the "calorimeter" information
    fftjet::Grid2d<Real> calo(nEtaCells, caloEtaMin, caloEtaMax, nPhiCells, 0.0);

    // The array which will hold noise information. Note that in our
    // model the noise is constant in E but not in Et. We will need to
    // reflect this by providing eta-dependent Et noise.
    std::vector<BgData> noiseData;
    noiseData.reserve(nEtaCells);
    for (unsigned i=0; i<nEtaCells; ++i)
    {
        const double eta = calo.etaBinCenter(i);
        noiseData.push_back(caloNoise/cosh(eta));
    }

    // The vector of preclusters
    std::vector<fftjet::Peak> preclusters;

    // The clustering tree (will be regrown every event)
    fftjet::PeakEtaPhiDistance treeDistance(etaToPhiBandwidthRatio);
    ClusteringTree theTree(&treeDistance);

    // The vector of reconstructed jets (we will refill it in every event)
    std::vector<RecoJet> recoJets;

    // Unclustered energy (vector and scalar) which we will refill
    // in every event. The "unclusScalar" will be the scalar sum of
    // the grid values not assigned to any jet. Due to the way
    // energy discretization is handled in this example, it will be
    // the scalar Et sum.
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

        // Regrow the clustering tree
        int status = clusteringSequencer.run(calo, &theTree);
        if (status)
        {
            std::cerr << "Event " << iev
                      << " : ERROR in ClusteringSequencer, status "
                      <<  status << std::endl;
            continue;
        }

        // Run the pattern recognition algorithm
        PatternRecoInfo patternInfo;
        status = patternRecoAlg.run(theTree, &preclusters, &patternInfo);
        if (status)
        {
            std::cerr << "Event " << iev
                      << " : ERROR in pattern recognition, status "
                      <<  status << std::endl;
            continue;
        }

        // Iterate the energy recombination stage
        unsigned nIterationsPerformed;
        double finalSquaredDistance;
        status = iterateRecombinationAlg(
            recoAlg, preclusters, calo,
            &noiseData[0], noiseData.size(), 1,
            jetCorrector, maxIterations, etaToPhiBandwidthRatio,
            relativePtBandwidth, convergedDistanceSquared,
            &recoJets, &unclustered, &unclusScalar,
            &nIterationsPerformed, &finalSquaredDistance);
        if (status)
        {
            std::cerr << "Event " << iev
                      << " : ERROR in reco algorithm iterations, status "
                      <<  status << std::endl;
            continue;
        }

        // In principle, we are done here. We now have the vector
        // of reconstructed jets and the unclustered energy. Run
        // some post-processing steps and print the results. Also,
        // "patternInfo", "nIterationsPerformed", and "finalSquaredDistance
        // variables contain some interesting information as well, but
        // they are not printed here.
        processResults(iev, etaToPhiBandwidthRatio, mcEvent,
                       recoJets, unclustered, unclusScalar);
    }

    // Clean up the objects which were created on the heap
    delete jetResponse;
    delete jetMemFcn;

    return 0;
}
