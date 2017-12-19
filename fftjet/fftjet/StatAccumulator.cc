#include <cfloat>
#include <cassert>
#include <cmath>

#include "fftjet/StatAccumulator.hh"

namespace fftjet {
    void StatAccumulator::reset()
    {
        sum_ = 0.0L;
        sumsq_ = 0.0L;
        min_ = DBL_MAX;
        max_ = -DBL_MAX;
        count_ = 0;
    }

    StatAccumulator::StatAccumulator()
    {
        reset();
    }

    void StatAccumulator::accumulate(const double value)
    {
        const long double ldv = value;
        sum_ += ldv;
        sumsq_ += ldv*ldv;
        ++count_;
        if (value < min_)
            min_ = value;
        if (value > max_)
            max_ = value;
    }

    void StatAccumulator::accumulate(const StatAccumulator& a)
    {
        sum_ += a.sum_;
        sumsq_ += a.sumsq_;
        count_ += a.count_;
        if (a.min_ < min_)
            min_ = a.min_;
        if (a.max_ > max_)
            max_ = a.max_;
    }

    double StatAccumulator::mean() const
    {
        assert(count_);
        return sum_/count_;
    }

    double StatAccumulator::stdev() const
    {
        assert(count_);
        if (count_ == 1)
            return 0.0;
        const long double var = (sumsq_ - sum_/count_*sum_)/(count_ - 1);
        if (var > 0.0L)
            return sqrt(static_cast<double>(var));
        else
            return 0.0;
    }
}
