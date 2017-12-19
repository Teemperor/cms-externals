//=========================================================================
// Matrix.hh
//
// A simple helper class for matrix manipulations. Depending on how much
// space is provided with the "Len" parameter, the data will be placed
// either on the stack or on the heap.
//
// C storage convention is used for internal data. In the element access
// operations, array bounds are not checked.
//
// Note that this class is slower and less complete than dedicated matrix
// classes based on expression templates (such as those in boost uBLAS
// or in Blitz++). Don't use it for any calculations in which speed is
// really important.
//
// I. Volobouev
// November 2008
//=========================================================================

#ifndef FFTJET_MATRIX_HH_
#define FFTJET_MATRIX_HH_

#include <iostream>

namespace fftjet {
    template<typename Numeric, unsigned Len=4>
    class Matrix
    {
    public:
        // Default constructor creates an unitialized matrix
        // which can be assigned from other matrix
        Matrix();

        // The following constructor creates an unitialized matrix
        // which can be assigned element-by-element or from another
        // matrix with the same dimensions
        Matrix(unsigned nrows, unsigned ncols);

        // The following constructor initializes the matrix as follows:
        //   initCode = 0     All elements are initialized to 0
        //   initCode = 1     Matrix must be square; diagonal elements
        //                    are initialized to 1
        Matrix(unsigned nrows, unsigned ncols, int initCode);

        // The following constructor initializes the matrix from the
        // given 1d array using C storage conventions
        Matrix(unsigned nrows, unsigned ncols, const Numeric* data);

        Matrix(const Matrix&);
        ~Matrix();

        Matrix& operator=(const Matrix&);

        inline unsigned nrows() const {return nrows_;}
        inline unsigned ncols() const {return ncols_;}
        inline Numeric* data() const {return data_;}

        // The following function resets the object to an unintialized state
        void reset();

        // The following function changes the object dimensions.
        // All data is lost in the process.
        void resize(unsigned nrows, unsigned ncols);

        // The following function sets all elements to 0
        void zeroOut();

        bool operator==(const Matrix&) const;
        bool operator!=(const Matrix&) const;

        Numeric* operator[](unsigned) const;

        Matrix operator*(const Matrix& r) const;
        Matrix operator*(Numeric r) const;
        Matrix operator/(Numeric r) const;
        Matrix operator+(const Matrix& r) const;
        Matrix operator-(const Matrix& r) const;
        Matrix operator+() const;
        Matrix operator-() const;

        Matrix& operator*=(Numeric r);
        Matrix& operator/=(Numeric r);
        Matrix& operator+=(const Matrix& r);
        Matrix& operator-=(const Matrix& r);

        Matrix T() const;
        Matrix symmetrize() const;

    private:
        Numeric local_[Len];
        Numeric* data_;
        unsigned nrows_;
        unsigned ncols_;
        unsigned len_;
    };
}

template<typename N, unsigned Len>
std::ostream& operator<<(std::ostream& os, const fftjet::Matrix<N, Len>& m);

#include "fftjet/Matrix.icc"

#endif // FFTJET_MATRIX_HH_
