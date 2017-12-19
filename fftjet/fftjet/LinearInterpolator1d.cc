#include <algorithm>
#include <cassert>
#include <cmath>

#include "fftjet/LinearInterpolator1d.hh"
#include "fftjet/cdfCellNumber.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/binaryIO.hh"

namespace fftjet {
    LinearInterpolator1d::LinearInterpolator1d(const double c)
        : cdf(0),
          xmin(0.0),
          flow(c),
          fhi(c),
          binwidth(1.0),
          npoints(2)
    {
        data = new double[2];
        data[0] = c;
        data[1] = c;
    }

    LinearInterpolator1d::LinearInterpolator1d(const LinearInterpolator1d& r)
        : Functor1<double, double>(),
          cdf(0),
          xmin(r.xmin),
          flow(r.flow),
          fhi(r.fhi),
          binwidth(r.binwidth),
          npoints(r.npoints)
    {
        data = new double[npoints];
        for (unsigned i=0; i<npoints; ++i)
            data[i] = r.data[i];
    }

    LinearInterpolator1d::LinearInterpolator1d(
        const std::vector<std::pair<double,double> >& v_in,
        const unsigned inpoints)
        : cdf(0),
          npoints(inpoints)
    {
        const unsigned len = v_in.size();
        assert(len > 1);
        assert(npoints > 1);
        std::vector<std::pair<double,double> > v(v_in);
        std::sort(v.begin(), v.end());
        binwidth = (v[len-1].first - v[0].first)/(npoints - 1);
        xmin = v[0].first - binwidth/2.0;
        flow = v[0].second;
        fhi = v[len-1].second;
        data = new double[npoints];
        data[0] = flow;
        data[npoints - 1] = fhi;
        unsigned ibelow = 0, iabove = 1;
        for (unsigned i=1; i<npoints-1; ++i)
        {
            const double x = xmin + (i + 0.5)*binwidth;
            while (v[iabove].first <= x)
            {
                ++ibelow;
                ++iabove;
            }
            if (v[ibelow].first == v[iabove].first)
                data[i] = (v[ibelow].second + v[iabove].second)/2.0;
            else
                data[i] = lin_interpolate_1d(
                    v[ibelow].first, v[iabove].first,
                    v[ibelow].second, v[iabove].second, x);
        }
    }

    LinearInterpolator1d& LinearInterpolator1d::operator=(
        const LinearInterpolator1d& r)
    {
        if (this == &r)
            return *this;
        Functor1<double, double>::operator=(r);
        xmin = r.xmin;
        flow = r.flow;
        fhi = r.fhi;
        binwidth = r.binwidth;
        if (npoints != r.npoints)
        {
            delete [] data;
            npoints = r.npoints;
            data = new double[npoints];
        }
        for (unsigned i=0; i<npoints; ++i)
            data[i] = r.data[i];
        delete [] cdf;
        cdf = 0;
        return *this;
    }

    double LinearInterpolator1d::operator()(const double& x) const
    {
        if (x < xmin)
            return flow;
        else
        {
            const unsigned ux = static_cast<unsigned>((x - xmin)/binwidth);
            if (ux >= npoints)
                return fhi;
            else
            {
                const double delta = (x - xmin)/binwidth - ux;
                if (ux == 0)
                {
                    if (delta <= 0.5)
                        return data[0];
                    else
                    {
                        const double w1 = delta - 0.5;
                        return data[0]*(1.0 - w1) + data[1]*w1;
                    }
                }
                else if (ux == npoints - 1)
                {
                    if (delta >= 0.5)
                        return data[ux];
                    else
                    {
                        const double wless = 0.5 - delta;
                        return data[ux-1]*wless + data[ux]*(1.0 - wless);
                    }
                }
                else
                {
                    if (delta <= 0.5)
                    {
                        const double wless = 0.5 - delta;
                        return data[ux-1]*wless + data[ux]*(1.0 - wless);
                    }
                    else
                    {
                        const double w1 = delta - 0.5;
                        return data[ux]*(1.0 - w1) + data[ux+1]*w1;
                    }
                }
            }
        }
    }

    bool LinearInterpolator1d::isNonNegative() const
    {
        if (flow < 0.0)
            return false;
        if (fhi < 0.0)
            return false;
        for (unsigned i=0; i<npoints; ++i)
            if (data[i] < 0.0)
                return false;
        return true;        
    }

    double LinearInterpolator1d::integral() const
    {
        long double sum = 0.0L;
        for (unsigned i=0; i<npoints; ++i)
            sum += data[i];
        return sum*binwidth;
    }

    void LinearInterpolator1d::normalize(const double value)
    {
        const double integ = integral();
        assert(integ != 0.0);
        const double scale = value/integ;
        for (unsigned i=0; i<npoints; ++i)
            data[i] *= scale;
    }

    double LinearInterpolator1d::singleCellRandom(
        const double xlo, const double bwx,
        const double z0, const double z1, const double rnd)
    {
        const double a = (z1 - z0)/2.0;
        const double c = rnd*(z0 + z1)/2.0;
        const double x = invertQuadraticCdf(a, z0, c);
        return x*bwx + xlo;
    }

    double LinearInterpolator1d::random(const double rnd) const
    {
        if (cdf == 0)
            const_cast<LinearInterpolator1d*>(this)->buildCdf();
        assert(rnd >= 0.0 && rnd <= 1.0);
        const unsigned ix = cdfCellNumber(rnd, cdf, npoints+1);

        // Figure out where is the bin for which we
        // have to generate the random number
        double xlo, bwx, z0, z1;
        if (ix == 0)
        {
            // Left edge
            xlo = xmin;
            bwx = binwidth/2.0;
            z0 = data[0];
            z1 = data[0];
        }
        else if (ix < npoints)
        {
            // Middle zone
            xlo = xmin + (ix-0.5)*binwidth;
            bwx = binwidth;
            z0 = data[ix-1];
            z1 = data[ix];
        }
        else if (ix == npoints)
        {
            // Right edge
            xlo = xmin + (ix-0.5)*binwidth;
            bwx = binwidth/2.0;
            z0 = data[npoints-1];
            z1 = data[npoints-1];
        }
        else
            assert(0);

        // Figure out the fractional x excess
        double dx = lin_interpolate_1d(
            ix ? cdf[ix-1] : 0.0, cdf[ix], 0, 1, rnd);

        // Generate a random location inside this bin
        return singleCellRandom(xlo, bwx, z0, z1, dx);
    }

    void LinearInterpolator1d::buildCdf()
    {
        assert(cdf == 0);
        assert(isNonNegative());
        cdf = new double[npoints+1];
        long double sum = data[0];
        cdf[0] = sum;
        for (unsigned i=1; i<npoints; ++i)
        {
            sum += (data[i] + data[i-1]);
            cdf[i] = sum;
        }
        sum += data[npoints-1];
        const double norm = sum;
        assert(norm > 0.0);
        for (unsigned i=0; i<npoints; ++i)
            cdf[i] /= norm;
        cdf[npoints] = 1.0;            
    }

    double LinearInterpolator1d::getCdf(const double x) const
    {
        if (x <= xmin)
            return 0.0;
        else
        {
            const unsigned ux = static_cast<unsigned>((x - xmin)/binwidth);
            if (ux >= npoints)
                return 1.0;
            else
            {
                double z0, z1, w0, w1, frac;
                if (cdf == 0)
                    const_cast<LinearInterpolator1d*>(this)->buildCdf();
                const double delta = (x - xmin)/binwidth - ux;
                if (ux == 0)
                {
                    if (delta <= 0.5)
                    {
                        // Left corner
                        z0 = data[0];
                        z1 = data[0];
                        w0 = 0.0;
                        w1 = cdf[0];
                        frac = delta*2.0;
                    }
                    else
                    {
                        z0 = data[0];
                        z1 = data[1];
                        w0 = cdf[0];
                        w1 = cdf[1];
                        frac = delta-0.5;
                    }
                }
                else if (ux == npoints - 1)
                {
                    if (delta >= 0.5)
                    {
                        // Right corner
                        z0 = data[ux];
                        z1 = data[ux];
                        w0 = cdf[ux];
                        w1 = 1.0;
                        frac = (delta-0.5)*2.0;
                    }
                    else
                    {
                        z0 = data[ux-1];
                        z1 = data[ux];
                        w0 = cdf[ux-1];
                        w1 = cdf[ux];
                        frac = delta+0.5;
                    }
                }
                else
                {
                    if (delta <= 0.5)
                    {
                        z0 = data[ux-1];
                        z1 = data[ux];
                        w0 = cdf[ux-1];
                        w1 = cdf[ux];
                        frac = delta+0.5;
                    }
                    else
                    {
                        z0 = data[ux];
                        z1 = data[ux+1];
                        w0 = cdf[ux];
                        w1 = cdf[ux+1];
                        frac = delta-0.5;
                    }
                }
                const double zsum((z0 + z1)/2.0);
                if (zsum == 0.0)
                    return w0;
                else
                    return w0 + (w1 - w0)/zsum*frac*(z0 + (z1 - z0)/2.0*frac);
            }
        }
    }

    bool LinearInterpolator1d::write(std::ostream& of) const
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

        write_pod(of, xmin);
        write_pod(of, flow);
        write_pod(of, fhi);
        write_pod(of, binwidth);
        write_pod(of, npoints);

        assert(data);
        write_pod_array(of, data, npoints);

        const long here = of.tellp();
        of.seekp(plen);
        u = here - objstart;
        write_pod(of, u);
        of.seekp(here);

        return !of.bad() && !of.fail();
    }

    LinearInterpolator1d* LinearInterpolator1d::read(std::istream& in)
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

        double xmin;
        double flow;
        double fhi;
        double binwidth;
        unsigned npoints;

        read_pod(in, &xmin);
        read_pod(in, &flow);
        read_pod(in, &fhi);
        read_pod(in, &binwidth);
        read_pod(in, &npoints);

        LinearInterpolator1d* i = 0;
        if (in.good())
        {
            double* data = new double[npoints];
            read_pod_array(in, data, npoints);
            if (!in.bad() && !in.fail())
                i = new LinearInterpolator1d(data, npoints, xmin,
                                             xmin + npoints*binwidth,
                                             flow, fhi);
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

    bool LinearInterpolator1d::operator==(const LinearInterpolator1d& r) const
    {
        if (!(xmin == r.xmin &&
              flow == r.flow &&
              fhi == r.fhi &&
              npoints == r.npoints))
            return false;

        if (binwidth != 0.0 || r.binwidth != 0.0)
            if (fabs(binwidth - r.binwidth)/(fabs(binwidth) + fabs(r.binwidth)) > 1.e-14)
                return false;

        for (unsigned i=0; i<npoints; ++i)
            if (data[i] != r.data[i])
                return false;

        return true;
    }
}
