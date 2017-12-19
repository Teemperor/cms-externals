#include <cassert>

#include "fftjet/InterpolatedKernel.hh"
#include "fftjet/binaryIO.hh"

namespace fftjet {
    void InterpolatedKernel::unscaledRectangle(KernelSupportRectangle *r) const
    {
        r->xmin = in.xMin();
        r->xmax = in.xMax();
        r->ymin = in.yMin();
        r->ymax = in.yMax();
    }

    void InterpolatedKernel::unscaledRandom(double r1, double r2,
                                            double* p1, double* p2) const
    {
        in.random(r1, r2, p1, p2);
    }

    bool InterpolatedKernel::operator==(const InterpolatedKernel& r) const
    {
        return xScaleFactor() == r.xScaleFactor() &&
            yScaleFactor() == r.yScaleFactor() &&
            scalePower() == r.scalePower() &&
            in == r.in;
    }

    bool InterpolatedKernel::write(std::ostream& of) const
    {
        const long plen = of.tellp();
        unsigned u = 0;
        write_pod(of, u);

        const long objstart = of.tellp();
        u = classId();
        write_pod(of, u);

        u = version();
        write_pod(of, u);

        // Write out the data from the base
        double d = xScaleFactor();
        write_pod(of, d);
        d = yScaleFactor();
        write_pod(of, d);
        int i = scalePower();
        write_pod(of, i);

        in.write(of);

        const long here = of.tellp();
        of.seekp(plen);
        u = here - objstart;
        write_pod(of, u);
        of.seekp(here);

        return !of.bad() && !of.fail();
    }

    InterpolatedKernel* InterpolatedKernel::read(std::istream& in)
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

        double xs;
        read_pod(in, &xs);
        double ys;
        read_pod(in, &ys);
        int ipow;
        read_pod(in, &ipow);

        LinearInterpolator2d* k = LinearInterpolator2d::read(in);
        if (k == 0)
            return 0;

        const long here = in.tellg();
        if (in.bad() || in.fail() || 
            len != static_cast<unsigned>(here - objstart))
        {
            delete k;
            return 0;
        }

        InterpolatedKernel* kernel = new InterpolatedKernel(xs, ys, ipow, *k);
        delete k;
        return kernel;
    }
}
