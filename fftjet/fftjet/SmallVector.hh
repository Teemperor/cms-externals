//=========================================================================
// SmallVector.hh
//
// A vector class which maintains its data on the stack until its size
// exceeds the template parameter, after which it migrates to the heap.
// This is not a full-blown implementation, just something to keep daugher
// pointers in various tree nodes mainly on the stack. The contained type
// must possess the default constructor and the assignment operator
// (or to be a simple built-in type).
//
// I. Volobouev
// June 2010
//=========================================================================

#ifndef FFTJET_SMALLVECTOR_HH_
#define FFTJET_SMALLVECTOR_HH_

#include <cassert>

namespace fftjet {
    template<typename T, unsigned Len>
    class SmallVector
    {
    public:
        // Default constructor
        SmallVector();

        // Constructor which makes a vector with the given size
        SmallVector(unsigned n);

        // The copy constructor
        SmallVector(const SmallVector&);

        // Converting constructor. It looks more general than the copy
        // constructor, but the actual copy constructor has to be created
        // anyway -- otherwise the compiler will generate an incorrect
        // default copy constructor.
        template<unsigned Len2>
        SmallVector(const SmallVector<T,Len2>&);

        // Destructor
        ~SmallVector();

        // Assignment operator
        SmallVector& operator=(const SmallVector&);

        // Converting assignment operator
        template <unsigned Len2>
        SmallVector& operator=(const SmallVector<T,Len2>&);

        // Comparison for equality
        template <unsigned Len2>
        bool operator==(const SmallVector<T,Len2>& r) const;

        template <unsigned Len2>
        inline bool operator!=(const SmallVector<T,Len2>& r) const
            {return !(*this == r);}

        // Subscripting
        inline T& operator[](const unsigned i) {return data_[i];}
        inline const T& operator[](const unsigned i) const {return data_[i];}

        inline T& at(const unsigned i)
            {assert(i < size_); return data_[i];}
        inline const T& at(const unsigned i) const
            {assert(i < size_); return data_[i];}

        // Other members with obvious meaning
        inline unsigned size() const {return size_;}
        inline unsigned capacity() const {return capacity_;}
        inline bool empty() const {return !size_;}

        void clear();
        void reserve(unsigned nElements);
        void resize(unsigned nElements);

        void push_back(const T& elem);
        void insert(unsigned position, const T& elem);
        void erase(unsigned position);

        // The following method returns the vector size if the element
        // is not found, otherwise it returns the element's position
        unsigned find(const T& elem) const;

    private:
        template <typename T2, unsigned Len2>
        friend class SmallVector;

        T* makeBuffer(unsigned sizeNeeded);

        T localData_[Len];
        T* data_;
        unsigned capacity_;
        unsigned size_;
    };
}

#include "fftjet/SmallVector.icc"

#endif // FFTJET_SMALLVECTOR_HH_
