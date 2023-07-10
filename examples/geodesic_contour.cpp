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
    std::vector<VertexIndex> sources;
    std::generate_n(std::back_insert_iterator(sources), 3, [&tm] {
        return VertexIndex(
            CGAL::get_default_random().uniform_int<int>(0, tm.num_vertices()));
    });
    hm.add_sources(sources);

    spdlog::info("sources: {}", sources);
    hm.estimate_geodesic_distances(dist_prop);

    IsoContour iso_contour(tm, dist_prop);


    return 0;
}