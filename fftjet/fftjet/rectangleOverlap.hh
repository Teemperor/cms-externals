//=========================================================================
// rectangleOverlap.hh
//
// Helper functions for the signal discretization process
//
// I. Volobouev
// January 2008
//=========================================================================

#ifndef FFTJET_RECTANGLEOVERLAP_HH_
#define FFTJET_RECTANGLEOVERLAP_HH_

#include <vector>
#include <utility>

namespace fftjet {
    // The "intervalOverlap" function returns the length of the overlap
    // of two intervals. If x1_min > x1_max or x2_min > x2_max then
    // the arguments will be swapped internally.
    double intervalOverlap(double x1_min, double x1_max,
                           double x2_min, double x2_max);

    // The "rectangleRectangleOverlap" function returns the area
    // of the overlap of two rectangles whose sides are parallel
    // to the coordinate axes.
    //
    // If x1_min > x1_max then arguments will be swapped internally.
    // Same is true for all other such pairs.
    //
    double rectangleRectangleOverlap(double x1_min, double y1_min,
                                     double x1_max, double y1_max,
                                     double x2_min, double y2_min,
                                     double x2_max, double y2_max);

    // The "rectanglePolygonOverlap" function returns the area
    // of the overlap of a rectangle whose sides are parallel
    // to the coordinate axes with a convex polygon. Polygon
    // is defined by its sequence of vertices. Each vertex is
    // represented by std::pair<double, double> where the first
    // element of the pair defines the x coordinate and the second
    // defines the y. The vertex sequence is represented by 
    // std::vector of such pairs.
    //
    // If xrect_min > xrect_max then arguments will be swapped
    // internally. Same is true for y.
    //
    // Don't use this function in those cases where the
    // "rectangleRectangleOverlap" function can be used.
    // "rectangleRectangleOverlap" runs significantly faster.
    //
    double rectanglePolygonOverlap(
        double xrect_min, double yrect_min,
        double xrect_max, double yrect_max,
        const std::vector<std::pair<double, double> >& polygon);
}

#endif // FFTJET_RECTANGLEOVERLAP_HH_
