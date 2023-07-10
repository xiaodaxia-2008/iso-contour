#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Surface_mesh.h>


namespace ic
{
namespace PMP = CGAL::Polygon_mesh_processing;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point_3>;
using VertexIndex = Mesh::Vertex_index;
using HalfedgeIndex = Mesh::Halfedge_index;
using EdgeIndex = Mesh::Edge_index;
using FaceIndex = Mesh::Face_index;
using FaceLocation = PMP::Face_location<Mesh, K::FT>;

class IntervalWithId : public CGAL::Interval_skip_list_interval<double>
{
public:
    HalfedgeIndex h;
    IntervalWithId(HalfedgeIndex h_, double inf_, double sup_, bool lb = true,
                   bool rb = true)
        : CGAL::Interval_skip_list_interval<double>(inf_, sup_, lb, rb), h(h_)
    {
    }
};

struct IsoContour
{
    const Mesh &tm;
    Mesh::Property_map<VertexIndex, double> scalar_prop;
    CGAL::Interval_skip_list<IntervalWithId> isl;

    /**
     * @note to convert FaceLocation to a 3D point, use
     * PMP::construct_point(floc, surface_mesh)
     *
     */
    typedef std::vector<FaceLocation> ContourSegment;
    typedef std::vector<ContourSegment> Contour;

    IsoContour(const Mesh &tm,
               const Mesh::Property_map<VertexIndex, double> &scalar_prop);
    IsoContour(const Mesh &tm, const std::string &vertex_property_name);

    Contour operator()(double iso_val);
};

} // namespace ic