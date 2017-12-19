#include "fftjet/Peak.hh"

namespace fftjet {
    Peak::Peak(const double eta, const double phi, const double mag,
               const double hessian[3],
               const double drsp, const double magsp,
               const double lifetime, const double scale,
               const double nearestDistance, const double memScale,
               const double recoScale, const double recoScaleRatio,
               const double clusterRadius, const double clusterSeparation,
               const int code, const int status,
               ScaleSpaceKernel* membershipFunction)
        : eta_(eta),
          phi_(phi),
          magnitude_(mag),
          speed_(drsp),
          magSpeed_(magsp),
          lifetime_(lifetime),
          splitTime_(-1.0),
          mergeTime_(-1.0),
          scale_(scale),
          nearestD_(nearestDistance),
          membershipFactor_(memScale),
          recoScale_(recoScale),
          recoScaleRatio_(recoScaleRatio),
          clusterRadius_(clusterRadius),
          clusterSeparation_(clusterSeparation),
          code_(code),
          status_(status),
          memFcn_(membershipFunction)
    {
        hessian_[0] = hessian[0];
        hessian_[1] = hessian[1];
        hessian_[2] = hessian[2];
    }

    void Peak::hessian(double hessianArray[3]) const
    {
        hessianArray[0] = hessian_[0];
        hessianArray[1] = hessian_[1];
        hessianArray[2] = hessian_[2];
    }

    void Peak::setHessian(const double hessianArray[3])
    {
        hessian_[0] = hessianArray[0];
        hessian_[1] = hessianArray[1];
        hessian_[2] = hessianArray[2];
    }

    double Peak::hessianDeterminant() const
    {
        return hessian_[0]*hessian_[2] - hessian_[1]*hessian_[1];
    }

    double Peak::scaledLaplacian(const double bwEta, const double bwPhi) const
    {
        // Note that both scale_*scale_*hessian_[0] and
        // scale_*scale_*hessian_[2] should simultaneously peak
        // at the correct scale. The addition of the scale factors here
        // simply makes the peak magnitudes comparable (of course,
        // only if the bandwidth ratio used for clustering correctly
        // reflects the bandwidth ratio in the data).
        return scale_*scale_*(hessian_[0]*bwEta*bwEta +
                              hessian_[2]*bwPhi*bwPhi);
    }

    bool Peak::operator==(const Peak& r) const
    {
        const bool hessians_equal = 
            hessian_[0] == r.hessian_[0] &&
            hessian_[1] == r.hessian_[1] &&
            hessian_[2] == r.hessian_[2];
        return hessians_equal &&
               eta_ == r.eta_ && 
               phi_ == r.phi_ && 
               magnitude_ == r.magnitude_ &&
               speed_ == r.speed_ &&
               magSpeed_ == r.magSpeed_ &&
               lifetime_ == r.lifetime_ &&
               splitTime_ == r.splitTime_ &&
               mergeTime_ == r.mergeTime_ &&
               scale_ == r.scale_ &&
               nearestD_ == r.nearestD_ &&
               membershipFactor_ == r.membershipFactor_ &&
               recoScale_ == r.recoScale_ &&
               recoScaleRatio_ == r.recoScaleRatio_ &&
               clusterRadius_ == r.clusterRadius_ &&
               clusterSeparation_ == r.clusterSeparation_ &&
               code_ == r.code_ &&
               status_ == r.status_;
        // There isn't any good way to compare membership functions
        // for equality. Comparing their pointers for identity does
        // not make sense here.
    }

    Peak Peak::dummy()
    {
        double hess[3] = {0.0, 0.0, 0.0};
        return Peak(0.0, 0.0, 0.0, hess);
    }
}
