//=========================================================================
// binaryIO.hh
//
// Utility functions for reading/writing "plain old data" in binary form
// (not platform independent yet)
//
// I. Volobouev
// April 2009
//=========================================================================

#ifndef FFTJET_BINARYIO_HH_
#define FFTJET_BINARYIO_HH_

#include <cassert>
#include <iostream>

namespace fftjet {
    // The following functions perform binary I/O of built-in types
    template <typename T>
    inline void write_pod(std::ostream& of, const T& pod)
    {
        of.write(reinterpret_cast<const char*>(&pod), sizeof(T));
    }
    
    template <typename T>
    inline void read_pod(std::istream& in, T* pod)
    {
        assert(pod);
        in.read(reinterpret_cast<char*>(pod), sizeof(T));
    }

    template <typename T>
    inline void write_pod_array(std::ostream& of, const T* pod,
                                const unsigned len)
    {
        if (len)
        {
            assert(pod);
            of.write(reinterpret_cast<const char*>(pod), len*sizeof(T));
        }
    }

    template <typename T>
    inline void read_pod_array(std::istream& in, T* pod, const unsigned len)
    {
        if (len)
        {
            assert(pod);
            in.read(reinterpret_cast<char*>(pod), len*sizeof(T));
        }
    }
}

#endif // FFTJET_BINARYIO_HH_
