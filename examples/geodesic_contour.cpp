#include "iso_contour.h"

#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <spdlog/fmt/ostr.h>
#include <spdlog/fmt/ranges.h>
#include <spdlog/spdlog.h>

using namespace ic;


template <>
struct fmt::formatter<VertexIndex> : ostream_formatter
{
};

int main()
{
    using HeatMethod =
        CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Mesh>;

    Mesh tm;
    const std::string fname = __FILE__ "/../tank.stl";
    if (!CGAL::IO::read_polygon_mesh(fname, tm)) {
        spdlog::error("failed to read polygon mesh from file {}", fname);
        return -1;
    }
    auto [dist_prop, created] =
        tm.add_property_map<VertexIndex, double>("v:distance", 0.0);

    HeatMethod hm(tm);
    hm.add_source(VertexIndex(2457));

    hm.estimate_geodesic_distances(dist_prop);

    spdlog::info("building interval tree...");
    IsoContour iso_contour(tm, dist_prop);
    spdlog::info("built interval tree");

    std::ofstream fvscalars(__FILE__ "/../vertex_scalar.txt");
    for (auto v : tm.vertices()) {
        fvscalars << fmt::format("{:.4f}\n", dist_prop[v]);
    }

    std::ofstream fcontours(__FILE__ "/../iso_contours.txt");
    for (auto contour_segment : iso_contour(0.3)) {
        for (auto loc : contour_segment) {
            auto p = PMP::construct_point(loc, tm);
            fcontours << fmt::format("{:.4f} {:.4f} {:4.f}\n", p[0], p[1],
                                     p[2]);
        }
    }

    return 0;
}