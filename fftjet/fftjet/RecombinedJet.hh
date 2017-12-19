//=========================================================================
// RecombinedJet.hh
//
// This class stores the jet information produced by the energy
// recombination algorithms. It includes the clustered quantity
// (e.g., a 4-vector) together with the original precluster from
// the pattern recognition stage, as well as several other
// quantities. This information can be retrieved using the
// following methods:
//
//   precluster()  -- Returns the precluster.
//
//   vec()         -- Returns the result of the recombination
//                    (normally, a 4-vector).
//
//   ncells()      -- The weighted number of energy discretization grid
//                    cells contributing to this jet. Depending on the
//                    algorithm settings (in particular, the "dataCutoff"
//                    parameter of KernelRecombinationAlg and similar
//                    classes), this number may or may not coincide with
//                    the jet area.
//
//   etSum()       -- The weighted sum of Et (pt, etc.) values
//                    contributing to this jet.
//
//   centroidEta() -- (Pseudo)rapidity of the weighted Et centroid.
//
//   centroidPhi() -- Azimuthal angle of the weighted Et centroid.
//
//   etaWidth()    -- Weighted eta RMS of the Et cells contributing
//                    to this jet.
//
//   phiWidth()    -- Weighted phi RMS of the Et cells contributing
//                    to this jet.
//
//   etaPhiCorr()  -- The weighted correlation coefficient between
//                    eta and phi.
//
//   fuzziness()   -- This quantity characterizes how far away the
//                    jet is from all other jets. In the "fuzzy" mode
//                    this is the fractional Et error due to the
//                    uncertainty in assigning the grid cells to this
//                    jet, calculated according to the generalized
//                    binomial distribution model. In the "crisp" mode
//                    this quantity does not have a well-defined meaning,
//                    but it is also a dimensionless number which
//                    becomes close to 0 for well-separated jets.
//
// All FFTJet algorithms produce this information. User-developed
// algorithms are encouraged to follow the same convention.
//
// In addition, there are two user-settable quantities with obvious
// suggested meaning which can be retrieved using methods "pileup()"
// and "uncertainty()". These quantities are not calculated by FFTJet
// algorithms. Another user-settable quantity, convergenceDistance(),
// could be used to indicate jet position and momentum convergence
// if the jet reconstruction is performed iteratively.
//
// Note that, in this package, "etaWidth" and "phiWidth" are defined
// with respect to the (eta, phi) centroid, so the width is the smallest
// possible. If you want the width with respect to the jet direction,
// add in quadrature the angular distance from the jet direction to the
// direction of the centroid.
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_RECOMBINEDJET_HH_
#define FFTJET_RECOMBINEDJET_HH_

#include "fftjet/Peak.hh"

namespace fftjet {
    template <typename VectorLike>
    class RecombinedJet
    {
    public:
        inline RecombinedJet(const Peak& peak, const VectorLike& p4,
                             double ncells, double etSum,
                             double centroidEta, double centroidPhi,
                             double etaWidth, double phiWidth,
                             double etaPhiCorr, double fuzziness,
                             double pileup=0.0, double uncertainty=0.0,
                             double convergenceDistance=-1.0)
            : peak_(peak), vec_(p4), ncells_(ncells),
              etSum_(etSum), centroidEta_(centroidEta),
              centroidPhi_(centroidPhi), etaWidth_(etaWidth),
              phiWidth_(phiWidth), etaPhiCorr_(etaPhiCorr),
              fuzziness_(fuzziness), pileup_(pileup),
              uncertainty_(uncertainty),
              convergenceDistance_(convergenceDistance) {}

        // The following functions simply return the arguments
        // provided in the class constructor
        inline const Peak& precluster() const {return peak_;}
        inline const VectorLike& vec() const {return vec_;}
        inline double ncells() const {return ncells_;}
        inline double etSum() const {return etSum_;}
        inline double centroidEta() const {return centroidEta_;}
        inline double centroidPhi() const {return centroidPhi_;}
        inline double etaWidth() const {return etaWidth_;}
        inline double phiWidth() const {return phiWidth_;}
        inline double etaPhiCorr() const {return etaPhiCorr_;}
        inline double fuzziness() const {return fuzziness_;}

        // It is possible to change the precluster
        inline void setPrecluster(const Peak& peak) {peak_ = peak;}

        // Modifiers for user-defined quantities
        inline void setPileup(const double d) {pileup_ = d;}
        inline void setUncertainty(const double d) {uncertainty_ = d;}
        inline void setConvergenceDistance(const double d)
            {convergenceDistance_ = d;}

        // Corresponding accessors
        inline double pileup() const {return pileup_;}
        inline double uncertainty() const {return uncertainty_;}
        inline double convergenceDistance() const
            {return convergenceDistance_;}

        // Some useful utilities
        bool operator==(const RecombinedJet& r) const;
        bool operator!=(const RecombinedJet& r) const;

        // We will use the 4-vector pt as the "magnitude"
        double magnitude() const;

        // Sorting by the magnitude
        bool operator<(const RecombinedJet& r) const;
        bool operator>(const RecombinedJet& r) const;

        // We need to provide functions which allow this class to serve
        // together with "Peak" in various templated code. The angular
        // functions are not provided in this manner: the user might want
        // to use either rapidity or pseudorapidity as "eta", and we can't
        // make this choice ahead of time. Nevertheless, the corresponding
        // peak functions are forwarded with a changed method name.
        inline double peakEta() const {return peak_.eta();}
        inline double peakPhi() const {return peak_.phi();}
        inline double peakMagnitude() const {return peak_.magnitude();}
        inline double driftSpeed() const {return peak_.driftSpeed();}
        inline double magSpeed() const {return peak_.magSpeed();}
        inline double lifetime() const {return peak_.lifetime();}
        inline double splitTime() const {return peak_.splitTime();}
        inline double mergeTime() const {return peak_.mergeTime();}
        inline double scale() const {return peak_.scale();}
        inline double nearestNeighborDistance() const
            {return peak_.nearestNeighborDistance();}
        inline double membershipFactor() const
            {return peak_.membershipFactor();}
        inline double recoScale() const
            {return peak_.recoScale();}
        inline double recoScaleRatio() const
            {return peak_.recoScaleRatio();}
        inline double clusterRadius() const
            {return peak_.clusterRadius();}
        inline double clusterSeparation() const
            {return peak_.clusterSeparation();}
        inline int code() const {return peak_.code();}
        inline int status() const {return peak_.status();}
        inline void hessian(double hessianArray[3]) const
            {peak_.hessian(hessianArray);}
        inline double hessianDeterminant() const
            {return peak_.hessianDeterminant();}
        inline double laplacian() const {return peak_.laplacian();}
        inline double scaledLaplacian(double bwEta, double bwPhi) const
            {return peak_.scaledLaplacian(bwEta, bwPhi);}
        inline ScaleSpaceKernel* membershipFunction() const
            {return peak_.membershipFunction();}

        inline void setPeakEtaPhi(const double eta, const double phi)
            {peak_.setEtaPhi(eta, phi);}
        inline void setPeakMagnitude(const double d) {peak_.setMagnitude(d);}
        inline void setDriftSpeed(const double d) {peak_.setDriftSpeed(d);}
        inline void setMagSpeed(const double d) {peak_.setMagSpeed(d);}
        inline void setLifetime(const double d) {peak_.setLifetime(d);}
        inline void setSplitTime(const double d) {peak_.setSplitTime(d);}
        inline void setMergeTime(const double d) {peak_.setMergeTime(d);}
        inline void setScale(const double d) {peak_.setScale(d);}
        inline void setNearestNeighborDistance(const double d)
            {peak_.setNearestNeighborDistance(d);}
        inline void setMembershipFactor(const double f)
            {peak_.setMembershipFactor(f);}
        inline void setRecoScale(const double s)
            {peak_.setRecoScale(s);}
        inline void setRecoScaleRatio(const double d)
            {peak_.setRecoScaleRatio(d);}
        inline void setClusterRadius(const double d)
            {peak_.setClusterRadius(d);}
        inline void setClusterSeparation(const double d)
            {peak_.setClusterSeparation(d);}
        inline void setCode(const int c) {peak_.setCode(c);}
        inline void setStatus(const int s) {peak_.setStatus(s);}
        inline void setHessian(const double hessian[3])
            {peak_.setHessian(hessian);}
        inline void setMembershipFunction(ScaleSpaceKernel* f)
            {peak_.setMembershipFunction(f);}

        // Sometimes there is a real need to create a dummy jet...
        static RecombinedJet dummy();

    private:
        RecombinedJet();

        Peak peak_;
        VectorLike vec_;
        double ncells_;
        double etSum_;
        double centroidEta_;
        double centroidPhi_;
        double etaWidth_;
        double phiWidth_;
        double etaPhiCorr_;
        double fuzziness_;
        double pileup_;
        double uncertainty_;
        double convergenceDistance_;
    };
}

#include "fftjet/RecombinedJet.icc"

#endif // FFTJET_RECOMBINEDJET_HH_
