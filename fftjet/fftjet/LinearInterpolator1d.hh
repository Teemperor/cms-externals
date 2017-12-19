//=========================================================================
// LinearInterpolator1d.hh
//
// This class can be used to interpolate histograms or similar regularly
// spaced data in one dimension
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_LINEARINTERPOLATOR1D_HH_
#define FFTJET_LINEARINTERPOLATOR1D_HH_

#include <vector>
#include <utility>
#include <iostream>

#include "fftjet/SimpleFunctors.hh"

namespace fftjet {
    class LinearInterpolator1d : public Functor1<double, double>
    {
    public:
        // Constructor from a regularly spaced data. "f_low"
        // is the function value below x_min, and "f_hi" is
        // the function value above "x_max".
        template <typename Real>
        LinearInterpolator1d(const Real* data, unsigned npoints,
                             double x_min, double x_max,
                             double f_low=0.0, double f_hi=0.0);

        // Constructor from a list of points, not necessarily regularly
        // spaced. The first member of the pair is the x coordinate, the
        // second is the tabulated function value. The input list will be
        // interpolated to "npoints" internal points linearly.
        LinearInterpolator1d(const std::vector<std::pair<double,double> >& v,
                             unsigned npoints);

        // Constructor which builds a function returning the given constant
        explicit LinearInterpolator1d(double c);

        // Copy constructor
        LinearInterpolator1d(const LinearInterpolator1d&);

        virtual ~LinearInterpolator1d();

        LinearInterpolator1d& operator=(const LinearInterpolator1d&);

        bool operator==(const LinearInterpolator1d& r) const;
        inline bool operator!=(const LinearInterpolator1d& r) const
            {return !(*this == r);}

        inline unsigned nx() const {return npoints;}
        inline double xMin() const {return xmin;}
        inline double xMax() const {return xmin + npoints*binwidth;}
        inline double fLow() const {return flow;}
        inline double fHigh() const {return fhi;}
        inline const double* getData() const {return data;}

        virtual double operator()(const double& x) const;

        // The following function tells whether the function
        // values are non-negative (and, therefore, the function
        // can be treated as a density)
        bool isNonNegative() const;

        // The following function returns the function integral
        // inside the interval
        double integral() const;

        // The following function returns the value of
        // the integral as if the interpolator is treated
        // as a density between xmin and xmax
        double getCdf(double x) const;

        // The following function can be used to generate
        // random numbers according to the function represented
        // by the interpolator, assuming that the interpolated
        // function can be used as a density
        double random(double rnd) const;

        // The following function can be used to normalize the
        // function integral inside the interval to the given value
        void normalize(double value);

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static LinearInterpolator1d* read(std::istream& in);

        // Helper function for generating random numbers
        // inside a single cell
        static double singleCellRandom(double xlo, double bwx,
                                       double z00, double z01,
                                       double rnd);
    private:
        LinearInterpolator1d();

        void buildCdf();

        double* data;
        double* cdf;
        double xmin;
        double flow;
        double fhi;
        double binwidth;
        unsigned npoints;

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 6;}
    };
}

#include "fftjet/LinearInterpolator1d.icc"

#endif // FFTJET_LINEARINTERPOLATOR1D_HH_
