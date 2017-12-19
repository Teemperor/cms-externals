//
// This example shows how to write clustering trees to disk in a manner
// suitable for subsequent visualization by OpenDX.
//
// Most of the code provided here is simply a subset of the code from
// "multiScaleExample.cc". Instead of running the pattern recognition
// and recombinatipon stages, this example simply writes the clustering
// trees out.
//
// When this program is run, it creates a subdirectory in the current
// directory named "ClusteringTrees", and a set of files in it named
// "clustree_NNN.dx", where NNN is the event number. If OpenDX is properly
// installed on your system, you can view the files as follows:
//
// 1) cd ClusteringTrees
//
// 2) Copy the "view_tree_sequence.net" program from the "opendx"
//    subdirectory of FFTJet into the current directory.
//
// 3) Run "dx -program view_tree_sequence.net". This will pop up the
//    OpenDX Visual Program Editor.
//
// 4) Click on the "Execute on Change" item in the "Execute" menu of
//    the Visual Program Editor. This will bring up the "Clustering Tree"
//    window.
//
// 5) Double-click on the "Sequencer" box inside the program window of the
//    Visual Program Editor. This will pop up the Sequence Control window.
//
// 6) Using the forward and backward buttons in the Sequence Control
//    window, navigate between events. You can interact with the images
//    inside the "Clustering Tree" window by pressing the left mouse
//    button and moving the cursor around (by default, this will invoke
//    the virtual trackball). You can change the interaction mode by
//    accessing the "Mode" menu inside the "Options" menu of the
//    "Clustering Tree" window. For more details about interacting with
//    the images, please see the OpenDX manual.
//
// I. Volobouev
// May 2009

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

// The following two includes are needed for "mkdir"
#include <sys/stat.h>
#include <sys/types.h>

// FFTJet headers we will need
#include "fftjet/Grid2d.hh"
#include "fftjet/ClusteringSequencer.hh"
#include "fftjet/DiscreteGauss2d.hh"
#include "fftjet/PeakSelectors.hh"
#include "fftjet/FrequencyKernelConvolver.hh"
#include "fftjet/EquidistantSequence.hh"
#include "fftjet/ProximityClusteringTree.hh"
#include "fftjet/PeakEtaPhiDistance.hh"
#include "fftjet/JetProperty.hh"
#include "fftjet/OpenDXPeakTree.hh"
#include "fftjet/ClusteringTreeSparsifier.hh"

// Headers from the "examples" directory
#include "fftjet_typedefs.hh"
#include "SimpleGenerator.hh"
#include "SimpleEvent.hh"
#include "discretizeMCEvent.hh"

int main(const int argc, char** argv)
{
    // Check if tree sparsification is requested
    const bool sparsify = argc == 2 && std::string(argv[1]) == "-sparse";

    // Some useful typedefs which simplify subsequent typing
    typedef fftjet::ProximityClusteringTree<fftjet::Peak,long> ClusteringTree;

    // Configure event topology and other job parameters

    // The number of events to process
    const unsigned long nevents = 100;

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
    const double maxScale = 1.0;

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

    // Subdirectory into which the files will be written.
    // There will be one file per event, named
    // subdir/clustree_NNN.dx, where NNN is the event number.
    const char *subdir = "ClusteringTrees";

    // Which features of the clustering tree we are going to visualize?
    // Each node in the tree will be displayed by OpenDX with a glyph
    // which will be positioned in the eta-phi-log(scale) space. Besides
    // the location, each glyph has, essentially, two characteristics
    // which can be used to visualize the data: size and color. In the
    // example given here, scale*scale*magnitude (which is proportional
    // to the jet energy at large scales) will be mapped into the glyph
    // size and scale-normalized Hessian (which might be potentially
    // useful as a blob detector) will be mapped into color.
    fftjet::ScaledMagnitude2<fftjet::Peak> glyphSize;
    fftjet::ScaledHessianDet<fftjet::Peak> glyphColor;

    /////////////////////////////////////////////////////////////////////////

    // Now, we will start constructing the relevant objects
    const double etaBandwidth = sqrt(etaToPhiBandwidthRatio);
    const double phiBandwidth = 1.0/etaBandwidth;
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

    /////////////////////////////////////////////////////////////////////////

    // Create the MC jet generator
    SimpleGenerator generator(njets, pythiaPartonCode, 
                              minPt, maxPt, minEta, maxEta, minPhi, maxPhi);

    // The structure which will hold the MC event info
    SimpleEvent mcEvent;

    // The structure which will hold the "calorimeter" information
    fftjet::Grid2d<Real> calo(nEtaCells,caloEtaMin,caloEtaMax,nPhiCells,0.0);

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

    // Processor for sparsifying the trees and the sparse tree
    fftjet::ClusteringTreeSparsifier<fftjet::Peak,long> sparsifier;
    fftjet::SparseClusteringTree<fftjet::Peak,long> sparseTree;

    // Create the directory for the output files
    // (read, write and execute permission by user)
    mkdir(subdir, S_IRWXU);

    // The following object will format the clustering trees
    // so that they can be later readable by OpenDX
    fftjet::OpenDXPeakTree<long,fftjet::AbsClusteringTree> formatter(
        &glyphSize, &glyphColor, caloEtaMax);

    // The following object will instead format sparse clustering trees
    fftjet::OpenDXPeakTree<long,fftjet::SparseClusteringTree> sparseFormatter(
        &glyphSize, &glyphColor, caloEtaMax);

    // Cycle over events
    for (unsigned long iev=0; iev<nevents; ++iev)
    {
        // Fill the MC event info
        generator.generate(&mcEvent);

        // Discretize the particle energies
        discretizeMCEvent(mcEvent, &calo, magneticFieldB,
                          caloRadius, caloNoise, caloThreshold);

        // Regrow the clustering tree
        clusteringSequencer.run(calo, &theTree);

        // In this example, we will insert the whole event
        // into the clustering tree
        clusteringSequencer.insertCompleteEvent(minScale*0.5, calo, &theTree);

        // Create the file name
        std::ostringstream filename;
        filename << subdir << '/' << "clustree_" << iev << ".dx";

        // Open the file
        std::ofstream file(filename.str().c_str());
        if (!file)
        {
            std::cerr << "Failed to open file \"" << filename.str()
                      << "\". Exiting." << std::endl;
            return 1;
        }

        if (sparsify)
        {
            // Fill the sparse tree
            sparsifier.sparsify(theTree, &sparseTree);

            // Write out the sparse tree
            sparseFormatter.setTree(sparseTree, 1, iev);
            file << sparseFormatter << std::endl;
        }
        else
        {
            // Just dump the tree. The file will be closed automatically
            // when the variable "file" goes out of scope.
            formatter.setTree(theTree, 1, iev);
            file << formatter << std::endl;
        }
    }

    return 0;
}
