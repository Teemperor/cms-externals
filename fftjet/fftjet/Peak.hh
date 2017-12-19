//=========================================================================
// Peak.hh
//
// Utility class which holds initial cluster (precluster) parameters.
// PeakFinder, ProximityClusteringTree, ClusteringSequencer, etc. use
// this class for its calculations.
//
// Objects of this class hold the following information:
//
//   eta()              -- Normally, rapidity or pseudorapidity of the
//                         precluster. The exact meaning of this quantity
//                         depends on the choices made by the user during
//                         the energy discretization step.
//
//   phi()              -- Azimuthal angle of the precluster.
//
//   magnitude()        -- The "magnitude" of the precluster: the height
//                         of the peak of the discretized energy
//                         distribution convoluted with the kernel.
//
//   scale()            -- Pattern recognition resolution scale at which
//                         this precluster was found.
//
//   hessian(array)     -- Returns the peak Hessian matrix (w.r.t. eta
//                         and phi variables) as a 1d array. Since this
//                         matrix is symmetric, it is enough to return
//                         three values. The order is h[0][0], h[0][1],
//                         h[1][1]. Note that both Hessian and Laplacian
//                         are provided by the peak finder only when
//                         it operates in the "subcell" resolution mode.
//
//   hessianDeterminant() -- Returns the determinant of the peak Hessian
//                           matrix.
//
//   laplacian()        -- Returns the peak Laplacian.
//
//   scaledLaplacian(bwEta, bwPhi) -- Returns the scale-normalized peak
//                                    Laplacian which can be used for blob
//                                    detection in the Gaussian scale space.
//
// The following quantities are defined if the preclusters were arranged
// in a clustering tree, and the corresponding calculations were performed
// by the tree code (by default, all these quantities are calculated):
//
//   driftSpeed()       -- The speed with which the precluster location
//                         moves in the scale space:
//                         d distance/d log(scale). Here, "distance"
//                         is defined in terms of the user-selected
//                         distance function.
//
//   magSpeed()         -- The speed with which the precluster magnitude
//                         changes in the scale space:
//                         d log(magnitude)/d log(scale)
//
//   lifetime()         -- The "lifetime" of the precluster in the scale
//                         space. It is computed as
//                         log(high_scale) - log(low_scale)
//                         where "low_scale" and "high_scale" define
//                         the range of resolution scales for which the
//                         precluster exists as a feature of the energy
//                         distribution. Typically, the lifetime is traced
//                         from the smallest scale in the clustering tree
//                         to the scale where the precluster becomes
//                         a part of a larger precluster. It may be useful
//                         to check the precluster lifetime in case
//                         the pattern recognition was performed with
//                         a kernel which tends to produce spurious
//                         modes. Such spurious preclusters simply
//                         disappear at smaller scales without leaving
//                         any descendants, which results in their low
//                         lifetimes.
//
//   nearestNeighborDistance() -- The distance to the nearest precluster
//                                at the same resolution scale.
//
// The following quantities are optional. They become really meaningful
// only if the complete event is inserted at the bottom level of the
// clustering tree, and this is a computationally expensive operation.
//
//   clusterRadius()     -- The distance from the precluster location
//                          to the farmost daughter
//
//   clusterSeparation() -- The distance from the precluster location
//                          to the nearest daughter of a different
//                          precluster
//
// The following quantities are user-settable. They are not calculated
// by FFTJet algorithms.
//
//   code()        --   The suggested use of these quantities is to
//   status()           indicate some feature and/or status of the
//                      user-developed pattern recognition algorithm.
//
// The following quantities may be used to steer the behavior of the
// energy recombination codes:
//
//   recoScale()   --   The initial recombination scale which can be
//                      set during pattern recognition stage and then
//                      passed to the jet membership function. If this
//                      scale is not set, the recombination stage will
//                      use the precluster resolution scale instead.
//
//   recoScaleRatio()    -- Ratio of eta to phi recombination scale.
//                          Can be set on per-cluster basis (in particular,
//                          this is useful if the typical cluster shape
//                          is changing depending on the cluster location
//                          in the detector).
//
//   membershipFactor()  -- Can be used to multiply the jet membership
//                          function by a precluster-dependent factor.
//                          It may be useful to set this factor, for example,
//                          in proportion to jet energy if the membership
//                          function is a normalized energy profile, etc.
//                          It can also be used for fast trimming of
//                          preclusters between pattern recognition and
//                          recombination stages -- just set the membership
//                          function factor to 0 to get rid of a precluster.
//
//  membershipFunction() -- User-settable jet membership function.
//                          By default, this function is not set, and
//                          all jets are recombined using the same
//                          jet membership function provided when
//                          the recombination algorithm object is created.
//                          However, in pursuit of ultimate resolution,
//                          the user may introduce precluster-specific
//                          assumptions (for example, about jet flavor).
//                          In this case it makes sense to use different
//                          membership functions for each jet. Note that
//                          all such function must have consistent
//                          dimensionality: either all of them should
//                          be derived from "AbsMembershipFunction" which
//                          includes an additional energy parameter or all
//                          of them should be derived from  "AbsKernel2d"
//                          which operates only in the eta-phi space.
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_PEAK_HH_
#define FFTJET_PEAK_HH_

#include <cmath>
#include <cassert>
#include <climits>

#include "fftjet/ScaleSpaceKernel.hh"

namespace fftjet {
    class Peak
    {
    public:
        // This object will not own the "membershipFunction" pointer
        Peak(double eta, double phi, double mag,
             const double hessian[3], double driftSpeed=-1.0,
             double magSpeed=-5.0, double lifetime=-1.0,
             double scale=-1.0, double nearestDistance=-1.0,
             double membershipFactor=1.0, double recoScale=0.0,
             double recoScaleRatio=0.0,
             double clusterRadius=-1.0, double clusterSeparation=-1.0,
             int code=INT_MIN, int status=-1,
             ScaleSpaceKernel* membershipFunction=0);

        // Inspectors
        inline double eta() const {return eta_;}
        inline double phi() const {return phi_;}
        inline double magnitude() const {return magnitude_;}
        inline double driftSpeed() const {return speed_;}
        inline double magSpeed() const {return magSpeed_;}
        inline double lifetime() const {return lifetime_;}
        inline double splitTime() const {return splitTime_;}
        inline double mergeTime() const {return mergeTime_;}
        inline double scale() const {return scale_;}
        inline double nearestNeighborDistance() const {return nearestD_;}
        inline double membershipFactor() const {return membershipFactor_;}
        inline double recoScale() const {return recoScale_;}
        inline double recoScaleRatio() const {return recoScaleRatio_;}
        inline double clusterRadius() const {return clusterRadius_;}
        inline double clusterSeparation() const {return clusterSeparation_;}
        inline int code() const {return code_;}
        inline int status() const {return status_;}
        void hessian(double hessianArray[3]) const;
        double hessianDeterminant() const;
        inline double laplacian() const {return hessian_[0] + hessian_[2];}
        double scaledLaplacian(double bwEta, double bwPhi) const;
        inline ScaleSpaceKernel* membershipFunction() const
            {return memFcn_;}

        // Modifiers
        inline void setEtaPhi(const double eta, const double phi)
            {eta_ = eta; phi_ = phi;}
        inline void setMagnitude(const double mag) {magnitude_ = mag;}
        inline void setDriftSpeed(const double s) {speed_ = s;}
        inline void setMagSpeed(const double s) {magSpeed_ = s;}
        inline void setLifetime(const double t) {lifetime_ = t;}
        inline void setSplitTime(const double t) {splitTime_ = t;}
        inline void setMergeTime(const double t) {mergeTime_ = t;}
        inline void setScale(const double scale) {scale_ = scale;}
        inline void setNearestNeighborDistance(const double d)
            {nearestD_ = d;}
        inline void setMembershipFactor(const double s)
            {membershipFactor_ = s;}
        inline void setRecoScale(const double d)
            {recoScale_ = d;}
        inline void setRecoScaleRatio(const double d)
            {recoScaleRatio_ = d;}
        inline void setClusterRadius(const double d)
            {clusterRadius_ = d;}
        inline void setClusterSeparation(const double d)
            {clusterSeparation_ = d;}
        inline void setCode(const int c) {code_= c;}
        inline void setStatus(const int s) {status_= s;}
        void setHessian(const double hessian[3]);
        inline void setMembershipFunction(ScaleSpaceKernel* f)
            {memFcn_ = f;}

        // Default comparison is by magnitude times scale squared
        inline bool operator<(const Peak& r) const
            {return scale_*scale_*magnitude_ < r.scale_*r.scale_*r.magnitude_;}

        inline bool operator>(const Peak& r) const
            {return scale_*scale_*magnitude_ > r.scale_*r.scale_*r.magnitude_;}

        bool operator==(const Peak& r) const;
        inline bool operator!=(const Peak& r) const {return !(*this == r);}

        // Normally, the default constructor of this object
        // should not be used: a peak without coordinates makes
        // no sense. The following function can help when
        // there is no way to avoid building a dummy object.
        static Peak dummy();

    private:
        Peak();

        double eta_;
        double phi_;
        double magnitude_;
        double speed_;
        double magSpeed_;
        double lifetime_;
        double splitTime_;
        double mergeTime_;
        double scale_;
        double nearestD_;
        double membershipFactor_;
        double recoScale_;
        double recoScaleRatio_;
        double clusterRadius_;
        double clusterSeparation_;
        double hessian_[3];
        int code_;
        int status_;
        ScaleSpaceKernel* memFcn_;
    };
}

#endif // FFTJET_PEAK_HH_
