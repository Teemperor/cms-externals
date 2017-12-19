//=========================================================================
// AbsClusterIntegrator.hh
//
// This class provides an interface definition to classes which
// estimate Pt (or energy) of the jets using just the initial
// preclusters provided by the clustering tree.
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef ABSCLUSTERINTEGRATOR_HH_
#define ABSCLUSTERINTEGRATOR_HH_

#include "fftjet/AbsClusteringTree.hh"

// A generic class for storing integration results.
// Exact interpretation of fields is up to the classes
// derived from "AbsClusterIntegrator". However, the following
// guidelines should be observed:
//
//   value    should be set to the integral value
//
//   error    should be set to the integral error or to a negative
//            number in case the error is not evaluated
//
//   status   negative status value means some kind of error occurred
//            during the calculation. In this case "value" is
//            meaningless. Interpretation of non-negative status
//            values is up to derived classes.
//
//   chisq,   use these values if the derived class performs
//   ndof     a least squares fit
//
class ClusterIntegral
{
public:
    inline ClusterIntegral(double value, double error, int status,
                           double chisq=0.0, unsigned ndof=0)
        : value_(value), error_(error), chisq_(chisq),
          ndof_(ndof), status_(status) {}

    inline double value() const {return value_;}
    inline double error() const {return error_;}
    inline int status() const {return status_;}
    inline double chisq() const {return chisq_;}
    inline unsigned ndof() const {return ndof_;}

    inline void setValue(const double v) {value_ = v;}
    inline void setError(const double e) {error_ = e;}
    inline void setStatus(const int s) {status_ = s;}
    inline void setChisq(const double c) {chisq_ = c;}
    inline void setNdof(const unsigned u) {ndof_ = u;}

    // Scaling by a constant factor
    inline ClusterIntegral& operator*=(const double& d)
    {
        value_ *= d;
        error_ *= d;
        chisq_ *= (d*d);
        return *this;
    }

private:
    ClusterIntegral();

    double value_;
    double error_;
    double chisq_;
    unsigned ndof_;
    int status_;
};

template<typename Cluster, typename LevelInfo>
struct AbsClusterIntegrator
{
    virtual ~AbsClusterIntegrator() {}

    virtual ClusterIntegral operator()(
        const fftjet::AbsClusteringTree<Cluster,LevelInfo>& tree,
        const fftjet::TreeNodeId& id) = 0;
};

#endif // ABSCLUSTERINTEGRATOR_HH_
