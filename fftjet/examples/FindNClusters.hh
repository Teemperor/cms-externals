//
// Example pattern recognition code which selects a level
// with the given number of jets from the clustering tree
//
// I. Volobouev
// May 2009

#ifndef FINDNCLUSTERS_HH_
#define FINDNCLUSTERS_HH_

#include "fftjet/AbsPatternRecognitionAlg.hh"
#include "fftjet/AbsClusteringTree.hh"
#include "fftjet/JetMagnitudeMapper2d.hh"
#include "fftjet/Peak.hh"

// The auxiliary data structure which provides additional
// information about the quality/significance of the pattern found.
//
// The members are as follows:
//
// jetCount  -- Actual number of jets found. If the configuration with
//              the requested number of jets can not be located, the
//              algorithm will attempt to find a configuration with
//              one more jet, and then with a smaller number of jets.
//
// minLevel,  -- The minimum and the maximum levels with the requested
// maxLevel      number of jets.
//
// foundLevel -- The level from which preclusters were taken.
// 
struct PatternRecoInfo
{
    inline PatternRecoInfo()
        : jetCount(0), minLevel(0), maxLevel(0), foundLevel(0) {}

    unsigned jetCount;
    unsigned minLevel;
    unsigned maxLevel;
    unsigned foundLevel;
};

class FindNClusters : public fftjet::AbsPatternRecognitionAlg<PatternRecoInfo>
{
public:
    // The constructor arguments are as follows:
    //   
    //   njets  -- We will search for configurations with this number of jets 
    //
    //   factor -- This is the factor for estimating the jet momentum
    //             from the peak magnitude. A crude Pt estimate can be
    //             obtaned from the product of factor*scale*scale*magnitude
    //             because the peak magnitude as a function of scale
    //             is essentially the Laplace transform of the jet
    //             energy profile. The "factor" is needed to properly
    //             normalize the transform, and it should normally be
    //             set to the total number of cells in the discretization
    //             grid divided by the eta range of the grid.
    //
    //   filename  -- The file which contains the peak magnitude response
    //                curve (as a function of log(scale) and original jet Pt).
    //
    // After constructing the object, the user should call the "isValid()"
    // method in order to verify that the object internals were successfully
    // initialized.
    // 
    FindNClusters(unsigned njets, double factor, const char* filename);
    ~FindNClusters();

    int run(const fftjet::AbsClusteringTree<fftjet::Peak,long>& tree,
            std::vector<fftjet::Peak>* peaks, PatternRecoInfo* info);

    inline bool isValid() const {return corrector && njets;}

private:
    typedef fftjet::JetMagnitudeMapper2d<fftjet::Peak> Corrector;

    FindNClusters();
    FindNClusters(const FindNClusters&);
    FindNClusters& operator=(const FindNClusters&);

    const double integFactor;
    Corrector* corrector;
    const unsigned njets;
};

#endif // FINDNCLUSTERS_HH_
