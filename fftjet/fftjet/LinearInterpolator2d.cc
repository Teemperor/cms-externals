#include <cassert>
#include <cmath>
#include <cfloat>

#include "fftjet/LinearInterpolator2d.hh"
#include "fftjet/cdfCellNumber.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/binaryIO.hh"

namespace fftjet {
    LinearInterpolator2d::LinearInterpolator2d(const double c)
        : cdf(0),
          nxpoints(2),
          nypoints(2),
          xmin(0.0),
          xbinwidth(1.0),
          ymin(0.0),
          ybinwidth(1.0),
          fOutside(c),
          useFOutside(true)
    {
        const unsigned npoints = 4;
        data = new double[npoints];
        for (unsigned i=0; i<npoints; ++i)
            data[i] = c;
    }

    LinearInterpolator2d::LinearInterpolator2d(const LinearInterpolator2d& r)
        : cdf(0),
          nxpoints(r.nxpoints),
          nypoints(r.nypoints),
          xmin(r.xmin),
          xbinwidth(r.xbinwidth),
          ymin(r.ymin),
          ybinwidth(r.ybinwidth),
          fOutside(r.fOutside),
          useFOutside(r.useFOutside)
    {
        const unsigned npoints = nxpoints*nypoints;
        data = new double[npoints];
        for (unsigned i=0; i<npoints; ++i)
            data[i] = r.data[i];
        if (r.cdf)
        {
            const unsigned ncdf = (nxpoints+1)*(nypoints+1);
            cdf = new double[ncdf];
            for (unsigned i=0; i<ncdf; ++i)
                cdf[i] = r.cdf[i];
        }
    }

    LinearInterpolator2d& LinearInterpolator2d::operator=(
        const LinearInterpolator2d& r)
    {
        if (this == &r)
            return *this;
        const bool reallocate = nxpoints*nypoints != r.nxpoints*r.nypoints;
        const bool reallocate_cdf = (nxpoints+1)*(nypoints+1) != 
                                    (r.nxpoints+1)*(r.nypoints+1);
        nxpoints = r.nxpoints;
        nypoints = r.nypoints;
        xmin = r.xmin;
        xbinwidth = r.xbinwidth;
        ymin = r.ymin;
        ybinwidth = r.ybinwidth;
        fOutside = r.fOutside;
        useFOutside = r.useFOutside;
        const unsigned npoints = nxpoints*nypoints;
        if (reallocate)
        {
            delete [] data;
            data = new double[npoints];
        }
        for (unsigned i=0; i<npoints; ++i)
            data[i] = r.data[i];

        if (r.cdf)
        {
            const unsigned ncdf = (nxpoints+1)*(nypoints+1);
            if (reallocate_cdf)
            {
                delete [] cdf;
                cdf = 0;
            }
            if (cdf == 0)
                cdf = new double[ncdf];
            for (unsigned i=0; i<ncdf; ++i)
                cdf[i] = r.cdf[i];
        }
        else
        {
            delete [] cdf;
            cdf = 0;
        }

        return *this;
    }

    inline double LinearInterpolator2d::cornerIntegral(
        const unsigned idx) const
    {
        double v;
        if (useFOutside)
            v = (3.0*fOutside + data[idx])/4.0;
        else
            v = data[idx];
        return v*xbinwidth*ybinwidth/4.0;
    }

    double LinearInterpolator2d::operator()(const double& x,
                                            const double& y) const
    {
        // Bin number
        int ix = x < xmin ? -1 : (int)floor((x - xmin)/xbinwidth);
        int iy = y < ymin ? -1 : (int)floor((y - ymin)/ybinwidth);

        // Point coordinates relative to the bin center in
        // units of the bin width
        double dx = (x - xmin)/xbinwidth - ix - 0.5;
        double dy = (y - ymin)/ybinwidth - iy - 0.5;

        if (ix < 0 || iy < 0 || 
            ix >= static_cast<int>(nxpoints) || 
            iy >= static_cast<int>(nypoints))
        {
            // Arguments outside the rectangle
            if (useFOutside)
                return fOutside;
            else
            {
                // Bring the point to the closest rectangle edge
                if (ix < 0)
                {
                    ix = 0;
                    dx = -0.5;
                }
                else if (ix >= static_cast<int>(nxpoints))
                {
                    ix = nxpoints - 1;
                    dx = 0.5;
                }
                if (iy < 0)
                {
                    iy = 0;
                    dy = -0.5;
                }
                else if (iy >= static_cast<int>(nypoints))
                {
                    iy = nypoints - 1;
                    dy = 0.5;
                }
            }
        }

        // Process all possible boundary conditions
        double z00, z10, z01, z11;
        if (ix == 0 && iy == 0 && dx <= 0.0 && dy <= 0.0)
        {
            // Lower left corner
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = fOutside;
                z01 = fOutside;
                z11 = data[0];
                dx  = 1.0 + dx*2.0;
                dy  = 1.0 + dy*2.0;
            }
            else
                return data[0];
        }
        else if (ix == 0 && iy == static_cast<int>(nypoints-1) && 
                 dx <= 0.0 && dy >= 0.0)
        {
            // Upper left corner
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = data[iy];
                z01 = fOutside;
                z11 = fOutside;
                dx  = 1.0 + dx*2.0;
                dy *= 2.0;
            }
            else
                return data[iy];
        }
        else if (ix == static_cast<int>(nxpoints-1) && iy == 0 && 
                 dx >= 0.0 && dy <= 0.0)
        {
            // Lower right corner
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = fOutside;
                z01 = data[ix*nypoints];
                z11 = fOutside;
                dx *= 2.0;
                dy  = 1.0 + dy*2.0;
            }
            else
                return data[ix*nypoints];
        }
        else if (ix == static_cast<int>(nxpoints-1) && 
                 iy == static_cast<int>(nypoints-1) && 
                 dx >= 0.0 && dy >= 0.0)
        {
            // Upper right corner
            if (useFOutside)
            {
                z00 = data[nxpoints*nypoints-1];
                z10 = fOutside;
                z01 = fOutside;
                z11 = fOutside;
                dx *= 2.0;
                dy *= 2.0;
            }
            else
                return data[nxpoints*nypoints-1];
        }
        else if (ix == 0 && dx <= 0.0)
        {
            // Left edge
            dx  = 1.0 + dx*2.0;
            if (dy < 0.0)
            {
                dy += 1.0;
                --iy;
            }
            z10 = data[iy];
            z11 = data[iy+1];
            if (useFOutside)
            {
                z00 = fOutside;
                z01 = fOutside;
            }
            else
            {
                z00 = z10;
                z01 = z11;
            }
        }
        else if (ix == static_cast<int>(nxpoints-1) && dx >= 0.0)
        {
            // Right edge
            dx *= 2.0;
            if (dy < 0.0)
            {
                dy += 1.0;
                --iy;
            }
            z00 = data[ix*nypoints+iy];
            z01 = data[ix*nypoints+iy+1];
            if (useFOutside)
            {
                z10 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z10 = z00;
                z11 = z01;
            }
        }
        else if (iy == 0 && dy <= 0.0)
        {
            // Bottom edge
            dy = 1.0 + dy*2.0;
            if (dx < 0.0)
            {
                dx += 1.0;
                --ix;
            }
            z01 = data[ix*nypoints];
            z11 = data[(ix+1)*nypoints];
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = fOutside;
            }
            else
            {
                z00 = z01;
                z10 = z11;
            }
        }
        else if (iy == static_cast<int>(nypoints-1) && dy >= 0.0)
        {
            // Top edge
            dy *= 2.0;
            if (dx < 0.0)
            {
                dx += 1.0;
                --ix;
            }
            z00 = data[ix*nypoints+iy];
            z10 = data[(ix+1)*nypoints+iy];
            if (useFOutside)
            {
                z01 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z01 = z00;
                z11 = z10;
            }
        }
        else
        {
            // We are not at the edge
            if (dx < 0.0)
            {
                dx += 1.0;
                --ix;
            }
            if (dy < 0.0)
            {
                dy += 1.0;
                --iy;
            }
            z00 = data[ix*nypoints+iy];
            z10 = data[(ix+1)*nypoints+iy];
            z01 = data[ix*nypoints+iy+1];
            z11 = data[(ix+1)*nypoints+iy+1];
        }

        return (dx - 1.0)*(dy - 1.0)*z00 + dx*z10 + 
            dy*z01 + dx*dy*(z11-z01-z10);
    }

    double LinearInterpolator2d::integral() const
    {
        long double sum = 0.0L;
        double z00, z10, z01, z11, bwx, bwy;

        // Process all boundary conditions

        // Process corners
        unsigned corners[4];
        corners[0] = 0;
        corners[1] = nypoints-1;
        corners[2] = (nxpoints-1)*nypoints;
        corners[3] = nxpoints*nypoints-1;
        for (unsigned i=0; i<4; ++i)
            sum += cornerIntegral(corners[i]);

        // Process left and right edges
        bwx = xbinwidth/2.0;
        bwy = ybinwidth;
        for (unsigned iy=0; iy<nypoints-1; ++iy)
        {
            z10 = data[iy];
            z11 = data[iy+1];
            if (useFOutside)
            {
                z00 = fOutside;
                z01 = fOutside;
            }
            else
            {
                z00 = z10;
                z01 = z11;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;

            const unsigned ix = nxpoints-1;
            z00 = data[ix*nypoints+iy];
            z01 = data[ix*nypoints+iy+1];
            if (useFOutside)
            {
                z10 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z10 = z00;
                z11 = z01;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
        }

        // Process bottom and top edges
        bwx = xbinwidth;
        bwy = ybinwidth/2.0;
        for (unsigned ix=0; ix<nxpoints-1; ++ix)
        {
            z01 = data[ix*nypoints];
            z11 = data[(ix+1)*nypoints];
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = fOutside;
            }
            else
            {
                z00 = z01;
                z10 = z11;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;

            const unsigned iy = nypoints-1;
            z00 = data[ix*nypoints+iy];
            z10 = data[(ix+1)*nypoints+iy];
            if (useFOutside)
            {
                z01 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z01 = z00;
                z11 = z10;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
        }

        // Process internal points
        bwx = xbinwidth;
        bwy = ybinwidth;
        for (unsigned ix=0; ix<nxpoints-1; ++ix)
            for (unsigned iy=0; iy<nypoints-1; ++iy)
                sum += (data[ix*nypoints+iy] + 
                        data[(ix+1)*nypoints+iy] + 
                        data[ix*nypoints+iy+1] + 
                        data[(ix+1)*nypoints+iy+1])/4.0*bwx*bwy;

        return sum;
    }

    void LinearInterpolator2d::normalize(const double value)
    {
        const double integ = integral();
        assert(integ != 0.0);
        const double scale = value/integ;
        const unsigned npoints = nxpoints*nypoints;
        for (unsigned i=0; i<npoints; ++i)
            data[i] *= scale;
    }

    bool LinearInterpolator2d::isNonNegative() const
    {
        if (useFOutside)
            if (fOutside < 0.0)
                return false;
        const unsigned npoints = nxpoints*nypoints;
        for (unsigned i=0; i<npoints; ++i)
            if (data[i] < 0.0)
                return false;
        return true;
    }

    LinearInterpolator2d::~LinearInterpolator2d()
    {
        delete [] cdf;
        delete [] data;
    }

    void LinearInterpolator2d::random(const double rnd1, const double rnd2,
                                      double *p1, double *p2) const
    {
        if (cdf == 0)
            const_cast<LinearInterpolator2d*>(this)->buildCdf();
        assert(rnd1 >= 0.0 && rnd1 <= 1.0);
        assert(rnd2 >= 0.0 && rnd2 <= 1.0);
        const unsigned nycdf = nypoints+1;

        // Figure out the x cell number
        const unsigned ix = cdfCellNumber(rnd1, cdf+nypoints,
                                          nxpoints+1, nycdf);

        // Figure out the y cell number. Map rnd2 so that we can find
        // the correct cell using the given range of cdf values.
        const double loval = ix ? cdf[(ix-1)*nycdf+nypoints] : 0.0;
        const double hival = cdf[ix*nycdf+nypoints];
        const double mapped = lin_interpolate_1d(0, 1, loval, hival, rnd2);
        const unsigned iy = cdfCellNumber(mapped, cdf+ix*nycdf, nycdf);

        // Figure out the fractional x excess
        const double dx = lin_interpolate_1d(loval, hival, 0, 1, rnd1);

        // Figure out the fractional y excess
        const unsigned idx = ix*nycdf+iy;
        const double hiy = cdf[idx];
        const double loy = idx ? cdf[idx-1] : 0.0;
        const double dy = lin_interpolate_1d(loy, hiy, 0, 1, mapped);

        // Figure out where is the bin for which we
        // have to generate the random number
        double xlo, ylo, bwx, bwy;
        if (ix == 0)
        {
            // Left edge
            xlo = xmin;
            bwx = xbinwidth/2.0;
        }
        else if (ix < nxpoints)
        {
            // Middle zone
            xlo = xmin + (ix-0.5)*xbinwidth;
            bwx = xbinwidth;
        }
        else if (ix == nxpoints)
        {
            // Right edge
            xlo = xmin + (ix-0.5)*xbinwidth;
            bwx = xbinwidth/2.0;
        }
        else
            assert(0);
        if (iy == 0)
        {
            // Bottom edge
            ylo = ymin;
            bwy = ybinwidth/2.0;
        }
        else if (iy < nypoints)
        {
            // Middle zone
            ylo = ymin + (iy-0.5)*ybinwidth;
            bwy = ybinwidth;
        }
        else if (iy == nypoints)
        {
            // Top edge
            ylo = ymin + (iy-0.5)*ybinwidth;
            bwy = ybinwidth/2.0;
        }
        else
            assert(0);

        // Now, figure out the function values
        // at the corners of this bin
        const double z00 = operator()(xlo, ylo);
        const double z01 = operator()(xlo, ylo+bwy);
        const double z10 = operator()(xlo+bwx, ylo);
        const double z11 = operator()(xlo+bwx, ylo+bwy);

        // Generate a random location inside this bin
        singleCellRandom(xlo, ylo, bwx, bwy, z00, z01, z10, z11,
                         dx, dy, p1, p2);
    }

    // The following cdf is really linear, not 2d
    void LinearInterpolator2d::buildCdf()
    {
        assert(cdf == 0);
        assert(isNonNegative());
        const unsigned ncdf = (nxpoints+1)*(nypoints+1);
        cdf = new double[ncdf];
        unsigned icdf = 0;

        long double sum = 0.0L;
        double z00, z10, z01, z11, bwx, bwy;

        // Left bottom corner
        sum += cornerIntegral(0);
        cdf[icdf++] = sum;

        // Left edge
        bwx = xbinwidth/2.0;
        bwy = ybinwidth;
        for (unsigned iy=0; iy<nypoints-1; ++iy)
        {
            z10 = data[iy];
            z11 = data[iy+1];
            if (useFOutside)
            {
                z00 = fOutside;
                z01 = fOutside;
            }
            else
            {
                z00 = z10;
                z01 = z11;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
            cdf[icdf++] = sum;
        }

        // Left top corner
        sum += cornerIntegral(nypoints-1);
        cdf[icdf++] = sum;

        // Central region
        bwx = xbinwidth;
        for (unsigned ix=0; ix<nxpoints-1; ++ix)
        {
            // Bottom edge
            bwy = ybinwidth/2.0;
            z01 = data[ix*nypoints];
            z11 = data[(ix+1)*nypoints];
            if (useFOutside)
            {
                z00 = fOutside;
                z10 = fOutside;
            }
            else
            {
                z00 = z01;
                z10 = z11;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
            cdf[icdf++] = sum;

            // Middle
            bwy = ybinwidth;
            for (unsigned iy=0; iy<nypoints-1; ++iy)
            {
                sum += (data[ix*nypoints+iy] + 
                        data[(ix+1)*nypoints+iy] + 
                        data[ix*nypoints+iy+1] + 
                        data[(ix+1)*nypoints+iy+1])/4.0*bwx*bwy;
                cdf[icdf++] = sum;
            }

            // Top edge
            bwy = ybinwidth/2.0;
            const unsigned iy = nypoints-1;
            z00 = data[ix*nypoints+iy];
            z10 = data[(ix+1)*nypoints+iy];
            if (useFOutside)
            {
                z01 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z01 = z00;
                z11 = z10;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
            cdf[icdf++] = sum;
        }

        // Right bottom corner
        sum += cornerIntegral((nxpoints-1)*nypoints);
        cdf[icdf++] = sum;

        // Right edge
        bwx = xbinwidth/2.0;
        bwy = ybinwidth;
        const unsigned ix = nxpoints-1;
        for (unsigned iy=0; iy<nypoints-1; ++iy)
        {
            z00 = data[ix*nypoints+iy];
            z01 = data[ix*nypoints+iy+1];
            if (useFOutside)
            {
                z10 = fOutside;
                z11 = fOutside;
            }
            else
            {
                z10 = z00;
                z11 = z01;
            }
            sum += (z10 + z11 + z00 + z01)/4.0*bwx*bwy;
            cdf[icdf++] = sum;
        }

        // Right top corner
        sum += cornerIntegral(nxpoints*nypoints-1);
        cdf[icdf++] = sum;
        assert(icdf == ncdf);

        // Renormalize
        const double dsum = static_cast<double>(sum);
        assert(dsum > 0.0);
        for (unsigned i=0; i<ncdf; ++i)
            cdf[i] /= dsum;
    }

    void LinearInterpolator2d::singleCellRandom(
        const double xlo, const double ylo,
        const double bwx, const double bwy,
        const double z00, const double z01,
        const double z10, const double z11,
        const double rnd1, const double rnd2,
        double* p1, double* p2)
    {
        const double integ = (z10 + z11 + z00 + z01)/4.0;

        // Solve a x^2 + b x = rnd1*integ, with appropriate a and b
        const double a = (z10 + z11 - z00 - z01)/4.0;
        const double b = (z00 + z01)/2.0;
        const double x = invertQuadraticCdf(a, b, rnd1*integ);

        // Do the same thing for y
        const double yinteg = 2.0*a*x + b;
        const double ya = (z01 - z00 + x*(z00 + z11 - z01 - z10))/2.0;
        const double yb = z00 - x*z00 + x*z10;
        const double y = invertQuadraticCdf(ya, yb, rnd2*yinteg);

        *p1 = xlo + x*bwx;
        *p2 = ylo + y*bwy;
    }

    bool LinearInterpolator2d::operator==(const LinearInterpolator2d& r) const
    {
        if (!(nxpoints == r.nxpoints &&
              nypoints == r.nypoints &&
              xmin == r.xmin &&
              fabs(xbinwidth - r.xbinwidth)/(fabs(xbinwidth) + fabs(r.xbinwidth)) < 1.e-14 &&
              ymin == r.ymin &&
              fabs(ybinwidth - r.ybinwidth)/(fabs(ybinwidth) + fabs(r.ybinwidth)) < 1.e-14 &&
              fOutside == r.fOutside &&
              useFOutside == r.useFOutside))
            return false;

        const unsigned npoints = nxpoints*nypoints;
        for (unsigned i=0; i<npoints; ++i)
            if (data[i] != r.data[i])
                return false;

        return true;
    }

    bool LinearInterpolator2d::write(std::ostream& of) const
    {
        const long plen = of.tellp();
        unsigned u = 0;
        write_pod(of, u);

        const long objstart = of.tellp();
        u = classId();
        write_pod(of, u);

        u = version();
        write_pod(of, u);

        // Write endianity id word
        u = 0x01020304U;
        write_pod(of, u);

        write_pod(of, nxpoints);
        write_pod(of, nypoints);
        write_pod(of, xmin);
        write_pod(of, xbinwidth);
        write_pod(of, ymin);
        write_pod(of, ybinwidth);
        write_pod(of, fOutside);
        write_pod(of, useFOutside);

        assert(data);
        write_pod_array(of, data, nxpoints*nypoints);

        const long here = of.tellp();
        of.seekp(plen);
        u = here - objstart;
        write_pod(of, u);
        of.seekp(here);

        return !of.bad() && !of.fail();
    }

    LinearInterpolator2d* LinearInterpolator2d::read(std::istream& in)
    {
        unsigned len;
        read_pod(in, &len);

        const long objstart = in.tellg();
        unsigned u;
        read_pod(in, &u);
        if (u != classId())
            return 0;

        read_pod(in, &u);
        if (u != version())
            return 0;

        read_pod(in, &u);
        if (u != 0x01020304U)
            return 0;

        unsigned nxpoints;
        unsigned nypoints;
        double xmin;
        double xbinwidth;
        double ymin;
        double ybinwidth;
        double fOutside;
        bool useFOutside;

        read_pod(in, &nxpoints);
        read_pod(in, &nypoints);
        read_pod(in, &xmin);
        read_pod(in, &xbinwidth);
        read_pod(in, &ymin);
        read_pod(in, &ybinwidth);
        read_pod(in, &fOutside);
        read_pod(in, &useFOutside);

        LinearInterpolator2d* i = 0;
        if (in.good())
        {
            double* data = new double[nxpoints*nypoints];
            read_pod_array(in, data, nxpoints*nypoints);
            if (!in.bad() && !in.fail())
                i = new LinearInterpolator2d(data, nxpoints, xmin,
                                             xmin + nxpoints*xbinwidth,
                                             nypoints, ymin,
                                             ymin + nypoints*ybinwidth,
                                             fOutside, useFOutside);
            delete [] data;
        }

        if (i)
        {
            // Check that the record length is correct
            const long here = in.tellg();
            if (len != static_cast<unsigned>(here - objstart))
            {
                delete i;
                i = 0;
            }
        }
        return i;
    }
}
