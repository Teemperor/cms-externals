//=========================================================================
// LinearInterpolator2d.hh
//
// This class can be used to interpolate histograms or similar regularly
// spaced data in two dimensions
//
// I. Volobouev
// August 2008
//=========================================================================

#ifndef FFTJET_LINEARINTERPOLATOR2D_HH_
#define FFTJET_LINEARINTERPOLATOR2D_HH_

#include <vector>
#include <utility>
#include <iostream>

namespace fftjet {
    class LinearInterpolator2d
    {
    public:
        // Constructor from a regularly spaced data.
        // Parameter "fOutside" can be used to specify the function
        // value in case the argument is outside the grid rectangle.
        // Parameter "useFOutside" tells whether, indeed, "fOutside"
        // should be used to do this or whether the function should
        // return the grid value at the closest rectangle edge.
        //
        // The "data" pointer can be NULL. In this case all internal
        // interpolator data will be set to 0. It is assumed that "setData"
        // method will be later used to provide the actual data.
        //
        template <typename Real>
        LinearInterpolator2d(const Real* data,
                             unsigned nxpoints, double xmin, double xmax,
                             unsigned nypoints, double ymin, double ymax,
                             double fOutside=0.0, bool useFOutside=true);

        // Constructor which builds a function returning the given constant
        explicit LinearInterpolator2d(double c);

        // Copy constructor
        LinearInterpolator2d(const LinearInterpolator2d&);

        ~LinearInterpolator2d();

        LinearInterpolator2d& operator=(const LinearInterpolator2d&);

        bool operator==(const LinearInterpolator2d& r) const;
        inline bool operator!=(const LinearInterpolator2d& r) const
            {return !(*this == r);}

        inline unsigned nx() const {return nxpoints;}
        inline unsigned ny() const {return nypoints;}
        inline double xMin() const {return xmin;}
        inline double xMax() const {return xmin + nxpoints*xbinwidth;}
        inline double yMin() const {return ymin;}
        inline double yMax() const {return ymin + nypoints*ybinwidth;}
        inline double outsideValue() const {return fOutside;}
        inline bool usingOutsideValue() const {return useFOutside;}
        inline const double* getData() const {return data;}

        double operator()(const double& x, const double& y) const;

        // The following function tells whether the function
        // values are non-negative (and, therefore, the function
        // can be treated as a density)
        bool isNonNegative() const;

        // The following function returns the function integral
        // inside the rectangle
        double integral() const;

        // The following function can be used to generate
        // random numbers according to the function represented
        // by the interpolator, assuming that the interpolated
        // function can be used as a density
        void random(double rnd1, double rnd2,
                    double *p1, double *p2) const;

        // The following function can be used to change the
        // interpolated data "on the fly". The number of points
        // in the input array should be at least as large as nx()*ny().
        template <typename Real>
        void setData(const Real* data);

        // The following function can be used to normalize the
        // function integral inside the rectangle to the given value
        void normalize(double value);

        // The "write" function returns "true" if the object
        // is successfully written out.
        bool write(std::ostream& of) const;

        // The following function reads the object in.
        // Returns NULL in case of a problem.
        static LinearInterpolator2d* read(std::istream& in);

        // Helper function for generating random numbers
        // inside a single cell
        static void singleCellRandom(double xlo, double ylo,
                                     double bwx, double bwy,
                                     double z00, double z01,
                                     double z10, double z11,
                                     double rnd1, double rnd2,
                                     double* p1, double* p2);
    private:
        LinearInterpolator2d();

        void buildCdf();
        double cornerIntegral(unsigned idx) const;

        double* data;
        double* cdf;
        unsigned nxpoints;
        unsigned nypoints;
        double xmin;
        double xbinwidth;
        double ymin;
        double ybinwidth;
        double fOutside;
        bool useFOutside;

        // Version number for the I/O
        static inline unsigned version() {return 1;}

        // Class id for the I/O
        static inline unsigned classId() {return 5;}
    };
}

#include "fftjet/LinearInterpolator2d.icc"

#endif // FFTJET_LINEARINTERPOLATOR2D_HH_
