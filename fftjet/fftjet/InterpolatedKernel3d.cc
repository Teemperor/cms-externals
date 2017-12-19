#include <cassert>
#include <algorithm>

#include "fftjet/InterpolatedKernel3d.hh"
#include "fftjet/interpolate.hh"
#include "fftjet/binaryIO.hh"

namespace fftjet {
    void InterpolatedKernel3d::supportRectangle(
        double, KernelSupportRectangle *r) const
    {
        r->xmin = xmin_;
        r->xmax = xmax_;
        r->ymin = ymin_;
        r->ymax = ymax_;
    }

    void InterpolatedKernel3d::setScaleRatio(double)
    {
    }

    double InterpolatedKernel3d::operator()(const double x, const double y,
                                            const double scale) const
    {
        assert(scale > 0.0);
        const unsigned nscales = scales_.size();
        const unsigned lb = std::lower_bound(scales_.begin(),
                                             scales_.end(), scale) -
                            scales_.begin();
        if (lb == 0)
            return interpols_[0](x, y);
        else if (lb == nscales)
            return interpols_[nscales-1](x, y);
        else
        {
            const unsigned lbm1(lb - 1);
            const double loval = interpols_[lbm1](x, y);
            const double hival = interpols_[lb](x, y);
            if (useLogSpace_)
                return log_interpolate_1d(scales_[lbm1], scales_[lb],
                                          loval, hival, scale);
            else
                return lin_interpolate_1d(scales_[lbm1], scales_[lb],
                                          loval, hival, scale);
        }
    }

    void InterpolatedKernel3d::random(const double r1, const double r2,
                                      const double scale,
                                      double* px, double* py) const
    {
        assert(scale > 0.0);
        if (scale == randomizerScale_)
            randomizer_.random(r1, r2, px, py);
        else
        {
            const unsigned nscales = scales_.size();
            const unsigned lb = std::lower_bound(scales_.begin(),
                                                 scales_.end(), scale) -
                scales_.begin();
            if (lb == 0)
                interpols_[0].random(r1, r2, px, py);
            else if (lb == nscales)
                interpols_[nscales-1].random(r1, r2, px, py);
            else
            {
                // Generating random numbers according to an arbitrary
                // distribution is very expensive. We need to cache some
                // intermediate results. First, check whether the scale
                // argument is close to some existing scale.
                const double eps = 1.e-8;
                const unsigned lbm1(lb - 1);
                const double loscale(scales_[lbm1]);
                const double hiscale(scales_[lb]);
                const double interval = hiscale - loscale;
                if (fabs(scale - loscale)/interval < eps)
                    interpols_[lbm1].random(r1, r2, px, py);
                else if (fabs(scale - hiscale)/interval < eps)
                    interpols_[lb].random(r1, r2, px, py);
                else
                {
                    buildRandomizer(lb, scale);
                    randomizerScale_ = scale;
                    randomizer_.random(r1, r2, px, py);
                }
            }
        }
    }

    void InterpolatedKernel3d::buildRandomizer(const unsigned lb,
                                               const double scale) const
    {
        const double loscale(scales_[lb - 1]);
        const double hiscale(scales_[lb]);
        const double whi = useLogSpace_ ?
            log(scale/loscale)/log(hiscale/loscale) :
            (scale - loscale)/(hiscale - loscale);
        const double wlo = 1.0 - whi;

        const double* zlo = interpols_[lb - 1].getData();
        const double* zhi = interpols_[lb].getData();
        const unsigned nbins = nxpoints_*nypoints_;
        double* newdata = new double[nbins];
        for (unsigned i=0; i<nbins; ++i)
            newdata[i] = wlo*zlo[i] + whi*zhi[i];
        randomizer_.setData(newdata);
        delete [] newdata;
    }

    bool InterpolatedKernel3d::isDensity() const
    {
        if (isDensity_ == 0)
        {
            isDensity_ = 1;
            const unsigned nscales = interpols_.size();
            for (unsigned i=0; i<nscales; ++i)
            {
                if (!interpols_[i].isNonNegative())
                {
                    isDensity_ = -1;
                    break;
                }
                if (interpols_[i].integral() <= 0.0)
                {
                    isDensity_ = -1;
                    break;
                }                
            }
        }
        return isDensity_ > 0;
    }

    bool InterpolatedKernel3d::operator==(const InterpolatedKernel3d& r) const
    {
        if (!(xmin_ == r.xmin_ &&
              xmax_ == r.xmax_ &&
              ymin_ == r.ymin_ &&
              ymax_ == r.ymax_ &&
              nxpoints_ == r.nxpoints_ &&
              nypoints_ == r.nypoints_ &&
              useLogSpace_ == r.useLogSpace_ &&
              scales_ == r.scales_))
            return false;

        if (!(interpols_ == r.interpols_))
            return false;

        // All other members are temporaries and don't affect the functionality
        return true;
    }

    bool InterpolatedKernel3d::write(std::ostream& of) const
    {
        // Write order:
        //
        // total length (not including this word)  -- unsigned
        // class id                                -- unsigned
        // version #                               -- unsigned
        // endianity identification word           -- unsigned
        //
        // everything else in a convenient order
        //
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

        write_pod(of, xmin_);
        write_pod(of, xmax_);
        write_pod(of, ymin_);
        write_pod(of, ymax_);
        write_pod(of, nxpoints_);
        write_pod(of, nypoints_);
        write_pod(of, useLogSpace_);

        u = scales_.size();
        write_pod(of, u);
        write_pod_array(of, &scales_[0], u);

        assert(u == interpols_.size());
        const unsigned npoints = nxpoints_*nypoints_;
        for (unsigned i=0; i<u; ++i)
            write_pod_array(of, interpols_[i].getData(), npoints);

        const long here = of.tellp();
        of.seekp(plen);
        u = here - objstart;
        write_pod(of, u);
        of.seekp(here);

        return !of.bad() && !of.fail();
    }

    InterpolatedKernel3d* InterpolatedKernel3d::read(std::istream& in)
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
        double xmax;
        double ymin;
        double ymax;
        unsigned nxpoints;
        unsigned nypoints;
        bool useLogSpace;

        read_pod(in, &xmin);
        read_pod(in, &xmax);
        read_pod(in, &ymin);
        read_pod(in, &ymax);
        read_pod(in, &nxpoints);
        read_pod(in, &nypoints);
        read_pod(in, &useLogSpace);

        read_pod(in, &u);
        std::vector<double> scales;
        scales.reserve(u);
        for (unsigned i=0; i<u; ++i)
        {
            double d;
            read_pod(in, &d);
            scales.push_back(d);
        }

        InterpolatedKernel3d* kernel = 0;
        if (in.good())
        {
            kernel = new InterpolatedKernel3d(
                static_cast<double*>(0), scales, useLogSpace, nxpoints,
                xmin, xmax, nypoints, ymin, ymax);

            const unsigned npoints = nxpoints*nypoints;
            double *data = new double[npoints];
            for (unsigned i=0; i<u; ++i)
            {
                read_pod_array(in, data, npoints);
                kernel->setScaleData(i, data);
            }
            delete [] data;
        }

        const long here = in.tellg();
        if (in.bad() || in.fail() ||
            len != static_cast<unsigned>(here - objstart))
        {
            delete kernel;
            kernel = 0;
        }
        return kernel;
    }
}
