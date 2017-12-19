//=========================================================================
// SimpleFunctors.hh
//
// Base classes for a variety of functor-based calculations
//
// I. Volobouev
// March 2009
//=========================================================================

#ifndef FFTJET_SIMPLEFUNCTORS_HH_
#define FFTJET_SIMPLEFUNCTORS_HH_

namespace fftjet {
    template <typename Result, typename Arg1>
    struct Functor1
    {
        typedef Arg1 argument_type;
        typedef Result result_type;

        virtual ~Functor1() {}
        virtual Result operator()(const Arg1&) const = 0;
    };

    template <typename Result, typename Arg1, typename Arg2>
    struct Functor2
    {
        typedef Arg1 first_argument_type;
        typedef Arg2 second_argument_type;
        typedef Result result_type;

        virtual ~Functor2() {}
        virtual Result operator()(const Arg1&, const Arg2&) const = 0;
    };

    // A simple functor which returns its argument
    template <typename Result>
    struct Same : public Functor1<Result, Result>
    {
        inline Result operator()(const Result& a) const {return a;}
    };

    // Adaptation for using single-argument functors with simple
    // cmath-style functions. Do not use this struct as a base class.
    template <typename Result, typename Arg1>
    struct FcnFunctor1 : public Functor1<Result, Arg1>
    {
        inline explicit FcnFunctor1(Result (*fcn)(Arg1)) : fcn_(fcn) {}

        inline Result operator()(const Arg1& a) const {return fcn_(a);}

    private:
        FcnFunctor1();
        Result (*fcn_)(Arg1);
    };

    // Adaptation for using two-argument functors with simple
    // cmath-style functions. Do not use this struct as a base class.
    template <typename Result, typename Arg1, typename Arg2>
    struct FcnFunctor2 : public Functor2<Result, Arg1, Arg2>
    {
        inline explicit FcnFunctor2(Result (*fcn)(Arg1, Arg2)) : fcn_(fcn) {}

        inline Result operator()(const Arg1& x, const Arg2& y) const
            {return fcn_(x, y);}

    private:
        FcnFunctor2();
        Result (*fcn_)(Arg1, Arg2);
    };
}

#endif // FFTJET_SIMPLEFUNCTORS_HH_
