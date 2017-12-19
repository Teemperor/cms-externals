#include <cassert>
#include <cmath>

#include "fftjet/cdfCellNumber.hh"
#include "fftjet/quartic_lib.hh"

namespace fftjet {
    unsigned cdfCellNumber(const double cdf, const double* cdfdata,
                           const unsigned ncells, const unsigned stride)
    {
        assert(cdfdata);
        assert(cdf >= 0.0 && cdf <= 1.0);

        if (ncells == 1)
            return 0;
        if (cdf == 0.0)
            for (unsigned i=0; i<ncells; ++i)
                if (cdfdata[stride*i] > 0.0)
                    return i;
        if (cdf == 1.0)
        {
            for (int i=ncells-2; i>=0; --i)
                if (cdfdata[stride*i] < 1.0)
                    return i+1;
            return 0;
        }
        unsigned imin = 0, imax = ncells-2;
        if (cdf < cdfdata[stride*imin])
            return imin;
        if (cdf >= cdfdata[stride*imax])
            return imax+1;
        while (imax - imin > 1)
        {
            const unsigned itry = (imin + imax)/2;
            if (cdf >= cdfdata[stride*itry])
                imin = itry;
            else
                imax = itry;
        }
        return imax;
    }

    double invertQuadraticCdf(const double a,
                              const double b,
                              const double c)
    {
        assert(c >= 0.0);
        unsigned nsols;
        double x0;
        if (a == 0.0)
        {
            assert(b > 0.0);
            nsols = 1;
            x0 = c/b;
        }
        else
        {
            double x1;
            nsols = fftjet::quadratic(b/a, -c/a, &x0, &x1);
            if (nsols == 2)
            {
                if (x0 > 0.0 && x0 < 1.0)
                    nsols = 1;
                else if (x1 > 0.0 && x1 < 1.0)
                {
                    x0 = x1;
                    nsols = 1;
                }
                else
                    nsols = 0;
            }
            if (nsols == 0)
            {
                if (c < fabs(a + b - c))
                    x0 = 0.0;
                else
                    x0 = 1.0;
                nsols = 1;
            }
        }
        assert(nsols == 1);
        return x0;
    }
}
