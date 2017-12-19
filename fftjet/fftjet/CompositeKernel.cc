#include <cfloat>
#include <cassert>

#include "fftjet/CompositeKernel.hh"
#include "fftjet/cdfCellNumber.hh"
#include "fftjet/interpolate.hh"

namespace fftjet {
    CompositeKernel::CompositeKernel(const bool takePointerOwnership)
        : destroyKernelSet_(takePointerOwnership)
    {
    }

    CompositeKernel::CompositeKernel(
        const std::vector<std::pair<double,AbsKernel2d*> >& kernSet,
        const bool takePointerOwnership)
        : std::vector<std::pair<double,AbsKernel2d*> >(kernSet),
          destroyKernelSet_(takePointerOwnership)
    {
    }

    CompositeKernel::~CompositeKernel()
    {
        if (destroyKernelSet_)
        {
            const unsigned n = size();
            for (unsigned i=0; i<n; ++i)
                delete (*this)[i].second;
        }
    }

    double CompositeKernel::operator()(const double x, const double y,
                                       const double scale) const
    {
        const unsigned n = size();
        double sum = 0.0;
        for (unsigned i=0; i<n; ++i)
        {
            const std::pair<double,AbsKernel2d*>* p = &(*this)[i];
            sum += p->first*(*p->second)(x, y, scale);
        }
        return sum;
    }

    void CompositeKernel::setScaleRatio(const double r)
    {
        const unsigned n = size();
        for (unsigned i=0; i<n; ++i)
        {
            std::pair<double,AbsKernel2d*>* p = &(*this)[i];
            p->second->setScaleRatio(r);
        }
    }

    double CompositeKernel::rectangleAverage(const double x,
                                             const double y,
                                             const double scale,
                                             const double dx,
                                             const double dy) const
    {
        const unsigned n = size();
        double sum = 0.0;
        for (unsigned i=0; i<n; ++i)
        {
            const std::pair<double,AbsKernel2d*>* p = &(*this)[i];
            sum += p->first*p->second->rectangleAverage(x, y, scale, dx, dy);
        }
        return sum;
    }

    void CompositeKernel::supportRectangle(const double scale,
                                           KernelSupportRectangle *r) const
    {
        const unsigned n = size();
        if (n == 0)
        {
            // What should be returned when the number
            // of components is 0? Not clear, but returning
            // 0 support seems to be reasonable...
            r->xmin = 0.0;
            r->xmax = 0.0;
            r->ymin = 0.0;
            r->ymax = 0.0;
            return;
        }

        r->xmin = DBL_MAX;
        r->xmax = -DBL_MAX;
        r->ymin = DBL_MAX;
        r->ymax = -DBL_MAX;
        KernelSupportRectangle tmp;
        for (unsigned i=0; i<n; ++i)
        {
            const std::pair<double,AbsKernel2d*>* p = &(*this)[i];
            if (p->first)
            {
                p->second->supportRectangle(scale, &tmp);
                if (tmp.xmin < r->xmin)
                    r->xmin = tmp.xmin;
                if (tmp.xmax > r->xmax)
                    r->xmax = tmp.xmax;
                if (tmp.ymin < r->ymin)
                    r->ymin = tmp.ymin;
                if (tmp.ymax > r->ymax)
                    r->ymax = tmp.ymax;
            }
        }
    }

    void CompositeKernel::random(const double r1, const double r2,
                                 const double scale,
                                 double* px, double* py) const
    {
        const unsigned n = size();
        assert(n);
        if (n != cdf_.size())
            cdf_.resize(n);
        double sum = 0.0;
        for (unsigned i=0; i<n; ++i)
        {
            const double weight((*this)[i].first);
            assert(weight >= 0.0);
            sum += weight;
            cdf_[i] = sum;
        }
        assert(sum > 0.0);
        for (unsigned i=0; i<n; ++i)
            cdf_[i] /= sum;
        const unsigned iker = cdfCellNumber(r1, &cdf_[0], n);
        const double loval = iker ? cdf_[iker-1] : 0.0;
        const double hival = cdf_[iker];
        const double mapped = lin_interpolate_1d(loval, hival, 0, 1, r1);
        ((*this)[iker].second)->random(mapped, r2, scale, px, py);
    }

    bool CompositeKernel::isDensity() const
    {
        const unsigned n = size();
        if (n == 0)
            return false;
        double sum = 0.0;
        for (unsigned i=0; i<n; ++i)
        {
            const double weight((*this)[i].first);
            if (weight < 0.0)
                return false;
            if (!(*this)[i].second->isDensity())
                return false;
            sum += weight;
        }
        if (sum <= 0.0)
            return false;
        return true;
    }
}
