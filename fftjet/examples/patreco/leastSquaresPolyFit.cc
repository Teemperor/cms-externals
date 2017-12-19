// This code is notably not thread-safe

#include <cmath>
#include <cassert>

#include "leastSquaresPolyFit.hh"

template<typename Numeric>
inline static void allocateStaticMemory(Numeric **memPtr, unsigned *len,
                                        const unsigned newlen)
{
    if (newlen > *len)
    {
        delete [] *memPtr;
        *memPtr = new Numeric[newlen];
        *len = newlen;
    }
}

// Declaration for the Lapack subroutines used
extern "C" {
    void dgels_(char *TRANS, int *M, int *N, int *NRHS, double *A, int *LDA,
                double *B, int *LDB, double *WORK, int *LWORK, int *INFO,
                int lenTRANS);
    void dtrtri_(char *UPLO, char *DIAG, int *N, double *A, int *LDA,
                 int *INFO, int lenUPLO, int lenDIAG);
}

double leastSquaresPolyFit(const unsigned polyDegree,
                           const double *x, const double *y,
                           const unsigned npoints, double *coeffs,
                           double *constTermError)
{
    // Space for the design matrix and other arrays
    static double *mem = 0;
    static unsigned lenMem = 0;

    assert(x);
    assert(y);
    assert(coeffs);
    assert(npoints > polyDegree);

    const unsigned degPlusOne = polyDegree + 1;

    // Allocate and parse out the memory
    int LWORK = (degPlusOne + 1)*32;
    allocateStaticMemory(&mem, &lenMem,
                         npoints*degPlusOne + npoints + LWORK);
    double *design = mem;
    double *values = design + npoints*degPlusOne;
    double *work   = values + npoints;

    // Build the design matrix
    if (polyDegree > 2)
        for (unsigned ideg=0; ideg+2<polyDegree; ++ideg)
        {
            const unsigned ipow = polyDegree - ideg;
            double *d = design + ideg*npoints;
            for (unsigned ipt=0; ipt<npoints; ++ipt)
                d[ipt] = pow(x[ipt], ipow);
        }
    if (polyDegree > 1)
    {
        double *d = design + (polyDegree-2)*npoints;
        for (unsigned ipt=0; ipt<npoints; ++ipt)
            d[ipt] = x[ipt]*x[ipt];
    }
    if (polyDegree > 0)
    {
        double *d = design + (polyDegree-1)*npoints;
        for (unsigned ipt=0; ipt<npoints; ++ipt)
            d[ipt] = x[ipt];
    }
    {
        double *d = design + polyDegree*npoints;
        for (unsigned ipt=0; ipt<npoints; ++ipt)
            d[ipt] = 1.0;
    }

    // Fill the vector of right hand sides
    for (unsigned ipt=0; ipt<npoints; ++ipt)
        values[ipt] = y[ipt];

    // Call the Lapack least squares routine
    char TRANS = 'N';
    int M = npoints;
    int N = degPlusOne;
    int NRHS = 1;
    int LDA = npoints;
    int LDB = npoints;
    int INFO = 1;
    dgels_(&TRANS, &M, &N, &NRHS, design, &LDA, values, &LDB,
           work, &LWORK, &INFO, 1);

    // Fill out the results
    assert(INFO == 0);
    for (unsigned ipow=0; ipow<degPlusOne; ++ipow)
        coeffs[ipow] = values[ipow];

    // Sum the squared residuals
    double dsumsq = 0.0;
    for (unsigned ipow=degPlusOne; ipow<npoints; ++ipow)
        dsumsq += values[ipow]*values[ipow];

    if (constTermError)
    {
        if (npoints > degPlusOne)
        {
            // Estimate of the constant term sigma.
            // Invert the upper triangular matrix.
            char UPLO = 'U';
            char DIAG = 'N';
            dtrtri_(&UPLO, &DIAG, &N, design, &LDA, &INFO, 1, 1);
            assert(INFO == 0);
            *constTermError = sqrt(dsumsq/(npoints - degPlusOne))*
                fabs(design[polyDegree*npoints+polyDegree]);
        }
        else
            *constTermError = 0.0;
    }

    return dsumsq;
}

double polyValue(const double *coeffs, unsigned polyDegree,
                 const double x)
{
    double res = *coeffs++;
    while (polyDegree--)
        res = res * x + *coeffs++;
    return res;
}
