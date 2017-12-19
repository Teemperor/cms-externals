//=========================================================================
// StatAccumulator.hh
//
// Simple single-pass (which means imprecise) accumulator of statistical
// summary for some set of numbers. Calculates minimum value, maximum
// value, mean, and standard deviation. "mean" and "stdev" functions
// will cause an assertion error in case "accumulate" function was never
// called after the object was created (or after it was reset).
//
// Use with care.
//
// I. Volobouev
// May 2008
//=========================================================================

#ifndef FFTJET_STATACCUMULATOR_HH_
#define FFTJET_STATACCUMULATOR_HH_

namespace fftjet {
    class StatAccumulator
    {
    public:
        StatAccumulator();

        inline unsigned long count() const {return count_;}
        inline double min() const {return min_;}
        inline double max() const {return max_;}        
        double mean() const;
        double stdev() const;

        void accumulate(double value);
        void accumulate(const StatAccumulator&);

        void reset();

    private:
        long double sum_;
        long double sumsq_;
        double min_;
        double max_;
        unsigned long count_;
    };
}

#endif // FFTJET_STATACCUMULATOR_HH_
