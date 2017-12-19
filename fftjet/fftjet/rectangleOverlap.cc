#include <cassert>
#include <cfloat>
#include <climits>
#include <algorithm>

#include "fftjet/rectangleOverlap.hh"
#include "fftjet/rectangleOverlap.scc"

namespace {
    // Simple helper class for the "rectanglePolygonOverlap" function.
    // This is a predicate for sorting points in the order of increasing
    // distance to the given point.
    class LesserDTo
    {
    public:
        inline LesserDTo(const std::pair<double, double>& input) : pt(input) {}
        inline bool operator()(const std::pair<double, double>& p1,
                               const std::pair<double, double>& p2) const
        {
            const double d1x(p1.first - pt.first);
            const double d1y(p1.second - pt.second);
            const double d2x(p2.first - pt.first);
            const double d2y(p2.second - pt.second);
            return d1x*d1x + d1y*d1y < d2x*d2x + d2y*d2y;
        }

    private:
        LesserDTo();
        const std::pair<double, double> pt;
    };
}

namespace fftjet {
    double intervalOverlap(double x1_min, double x1_max,
                           double x2_min, double x2_max)
    {
        swap_if_larger(x1_min, x1_max);
        swap_if_larger(x2_min, x2_max);

        double newmax = 0.0, newmin = 0.0;
        if (x1_max > x2_min && x2_max > x1_min)
        {
            newmin = x1_min < x2_min ? x2_min : x1_min;
            newmax = x1_max < x2_max ? x1_max : x2_max;
        }
        return newmax - newmin;
    }

    double rectangleRectangleOverlap(const double x1_min, const double y1_min,
                                     const double x1_max, const double y1_max,
                                     const double x2_min, const double y2_min,
                                     const double x2_max, const double y2_max)
    {
        double result = 0.0;
        const double ox(intervalOverlap(x1_min, x1_max, x2_min, x2_max));
        if (ox)
            result = ox*intervalOverlap(y1_min, y1_max, y2_min, y2_max);
        return result;
    }

    double rectanglePolygonOverlap(
        double xrect_min, double yrect_min,
        double xrect_max, double yrect_max,
        const std::vector<std::pair<double, double> >& polygon)
    {
        const unsigned n_vertices(polygon.size());
        if (n_vertices < 3)
            return 0.0;

        swap_if_larger(xrect_min, xrect_max);
        swap_if_larger(yrect_min, yrect_max);

        // First, check whether the input rectangle overlaps with
        // the rectangle that encloses the input polygon. This check
        // can be performed very quickly. If there is no overlap
        // then we are done.
        unsigned nv_inside_rectangle = 0;
        unsigned ivert_outside = UINT_MAX;
        {
            double xmin_poly = DBL_MAX, xmax_poly = -DBL_MAX;
            double ymin_poly = DBL_MAX, ymax_poly = -DBL_MAX;
            for (unsigned i=0; i<n_vertices; ++i)
            {
                const std::pair<double, double>& vertex(polygon[i]);
                if (vertex.first < xmin_poly)
                    xmin_poly = vertex.first;
                if (vertex.first > xmax_poly)
                    xmax_poly = vertex.first;
                if (vertex.second < ymin_poly)
                    ymin_poly = vertex.second;
                if (vertex.second > ymax_poly)
                    ymax_poly = vertex.second;
                if (vertex.first >= xrect_min && vertex.first <= xrect_max &&
                    vertex.second >= yrect_min && vertex.second <= yrect_max)
                    ++nv_inside_rectangle;
                else if (ivert_outside == UINT_MAX)
                    ivert_outside = i;
            }
            if (rectangleRectangleOverlap(xrect_min, yrect_min,
                                          xrect_max, yrect_max,
                                          xmin_poly, ymin_poly,
                                          xmax_poly, ymax_poly) == 0.0)
                return 0.0;

            // Make sure the polygon is convex
            assert(is_poly_convex(polygon));

            // Check whether the polygon is completely enclosed
            // by the rectangle
            if (nv_inside_rectangle == n_vertices)
                return fabs(polygon_area(polygon));
        }

        // No more useful simple tests.
        // Time to do the heavy lifting now.
        std::vector<std::pair<double, double> > rectangle;
        rectangle.reserve(4);
        rectangle.push_back(std::make_pair(xrect_min, yrect_min));
        rectangle.push_back(std::make_pair(xrect_min, yrect_max));
        rectangle.push_back(std::make_pair(xrect_max, yrect_max));
        rectangle.push_back(std::make_pair(xrect_max, yrect_min));

        // Check whether the rectangle is completely enclosed by the polygon
        unsigned n_rect_v_enclosed = 0;
        bool rect_vertice_enclosed[4];
        for (unsigned irect = 0; irect < 4; ++irect)
        {
            rect_vertice_enclosed[irect] = is_point_inside_poly(
                rectangle[irect].first, rectangle[irect].second, polygon);
            if (rect_vertice_enclosed[irect])
                ++n_rect_v_enclosed;
        }
        if (n_rect_v_enclosed == 4)
            return (xrect_max - xrect_min)*(yrect_max - yrect_min);

        // Find all the crossings of the rectangle and the polygon
        unsigned n_crossings = 0;
        unsigned first_crossing_at = 10;
        std::vector<std::vector<std::pair<double, double> > > crossings(4);
        std::vector<std::pair<double, double> > common;
        for (unsigned irect = 0; irect < 4; ++irect)
        {
            for (unsigned ipoly = 0; ipoly < n_vertices; ++ipoly)
                append_intersections(rectangle[irect], rectangle[(irect + 1) % 4],
                                     polygon[ipoly], polygon[(ipoly + 1) % n_vertices],
                                     &crossings[irect]);
            const unsigned n(crossings[irect].size());
            n_crossings += n;
            // Sort the crossing points so that they
            // will appear clockwise eventually
            if (n > 1)
                sort(crossings[irect].begin(), crossings[irect].end(),
                     LesserDTo(rectangle[irect]));
            if (n > 0 && first_crossing_at == 10)
            {
                first_crossing_at = irect;
                for (unsigned i=0; i<n; ++i)
                    common.push_back(crossings[irect][i]);
            }
        }
        if (n_crossings == 0)
        {
            assert(nv_inside_rectangle == 0);
            assert(n_rect_v_enclosed == 0);
            return 0.0;
        }

        assert(n_crossings > 1);
        assert(ivert_outside < n_vertices);

        // We will need to know whether the polygon was specified
        // in the clockwise or counterclockwise order. The simplest
        // (although not the fastest) way to do it here is to lookup
        // the sign of the polygon area. Note that the polygon area
        // at this point should not be 0 because convex polygons with
        // 0 area should be taken care of already by the
        // "rectangleRectangleOverlap" check at the beginning.
        const int area_s = sgn(polygon_area(polygon));
        assert(area_s);

        for (unsigned itmp = first_crossing_at+1; itmp < first_crossing_at+5; ++itmp)
        {
            const unsigned irect = itmp % 4;
            const unsigned n(crossings[irect].size());
            if (rect_vertice_enclosed[irect])
                common.push_back(rectangle[irect]);
            else if (n > 0 && nv_inside_rectangle > 0)
            {
                // We need to add the right polygon vertices
                // to the sequence in the clockwise order
                const unsigned ilast = common.size() - 1;
                if (area_s > 0)
                    for (unsigned i=ivert_outside+1; i<ivert_outside+n_vertices; ++i)
                    {
                        const std::pair<double, double>& vertex(polygon[i % n_vertices]);
                        if (vertex.first >= xrect_min && vertex.first <= xrect_max &&
                            vertex.second >= yrect_min && vertex.second <= yrect_max)
                            if (point_sign(vertex.first, vertex.second,
                                           common[ilast], crossings[irect][0]) < 0)
                                common.push_back(vertex);
                    }
                else
                    for (unsigned i=ivert_outside+n_vertices-1; i>ivert_outside; --i)
                    {
                        const std::pair<double, double>& vertex(polygon[i % n_vertices]);
                        if (vertex.first >= xrect_min && vertex.first <= xrect_max &&
                            vertex.second >= yrect_min && vertex.second <= yrect_max)
                            if (point_sign(vertex.first, vertex.second,
                                           common[ilast], crossings[irect][0]) < 0)
                                common.push_back(vertex);
                    }
            }
            if (itmp < first_crossing_at+4)
                for (unsigned i=0; i<n; ++i)
                    common.push_back(crossings[irect][i]);
        }

        return fabs(polygon_area(common));
    }
}
