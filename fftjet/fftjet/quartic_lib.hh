//=========================================================================
// quartic_lib.hh
//
// Numerically sound routines for solving algebraic equations.
// Original code by Don Herbison-Evans, with minimal adaptations for
// this package by igv. See the following webpage for a description
// of methods used:
//
// http://linus.socs.uts.edu.au/~don/pubs/solving.html
//
// I. Volobouev
// November 2008
//=========================================================================

#ifndef FFTJET_QUARTIC_LIB_HH_
#define FFTJET_QUARTIC_LIB_HH_

namespace fftjet {
    // Numerically stable code to solve quadratic equations:
    //   x**2 + b*x + c == 0
    // Returns the number of real roots found. The roots are placed 
    // into *x1 and *x2. The case of 0 discriminant is treated as if
    // there are two equal roots.
    int quadratic(double b, double c, double *x1, double *x2);

    // Find the real roots of the cubic:
    //   x**3 + p*x**2 + q*x + r = 0
    // The number of real roots is returned, and the roots are placed
    // into the "v3" array.
    int cubic(double p, double q, double r, double v3[3]);

    // Solve quartic equation using, hopefully, a numerically stable
    // scheme (see the web page for details). The number of real roots
    // is returned, and the roots are placed into the "rts" array.
    int quartic(double a, double b, double c, double d, double rts[4]);

    // Set the debug level for cubic and quartic.
    // <1 for lots of diagnostics, >5 for none.
    void set_quartic_debug(int level);
}

#endif /* FFTJET_QUARTIC_LIB_HH_ */
