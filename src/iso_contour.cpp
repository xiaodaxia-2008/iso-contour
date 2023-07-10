#include "iso_contour.h"

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <spdlog/spdlog.h>

namespace ic
{

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              FaceLocation>
    AL_graph;

struct Polyline_visitor
{
    std::vector<std::vector<FaceLocation>> polylines;
    AL_graph &al_graph;

    Polyline_visitor(AL_graph &al_graph_) : al_graph(al_graph_) {}

    void start_new_polyline() { polylines.push_back({}); }

    void add_node(AL_graph::vertex_descriptor node_id)
    {
        polylines.back().push_back(al_graph[node_id]);
    }

    void end_polyline()
    {
        if (polylines.back().empty()) {
            polylines.pop_back();
        }
    }
};

IsoContour::IsoContour(
    const Mesh &tm_,
    const Mesh::Property_map<VertexIndex, double> &scalar_prop_)
    : tm(tm_), scalar_prop(scalar_prop_)
{
    for (auto e : tm.edges()) {
        auto h = tm.halfedge(e);
        auto infimum = scalar_prop[tm.source(h)];
        auto supremum = scalar_prop[tm.target(h)];
        if (infimum > supremum) {
            std::swap(infimum, supremum);
            h = tm.opposite(h);
        }
        isl.insert(IntervalWithId(h, infimum, supremum));
    }
}

IsoContour::IsoContour(const Mesh &tm, const std::string &vertex_property_name)
    : IsoContour(
        tm, tm.property_map<VertexIndex, double>(vertex_property_name).first)
{
}

auto IsoContour::operator()(double iso_val) -> Contour
{
    std::vector<IntervalWithId> intervals;
    isl.find_intervals(iso_val, std::back_insert_iterator(intervals));
    SPDLOG_DEBUG("found {} intervals contain iso value {}", intervals.size(),
                 iso_val);

    auto CompareFace = [this](const HalfedgeIndex &lhs,
                              const HalfedgeIndex &rhs) {
        return this->tm.face(lhs) < this->tm.face(rhs);
    };
    typedef std::pair<AL_graph::vertex_descriptor, AL_graph::vertex_descriptor>
        AL_vertex_pair;
    std::map<HalfedgeIndex, AL_vertex_pair, decltype(CompareFace)> al_edge_map(
        CompareFace);

    AL_graph al_graph;
    for (const auto &interval : intervals) {
        auto h = interval.h;
        auto s0 = scalar_prop[tm.source(h)];
        auto s1 = scalar_prop[tm.target(h)];
        auto t = (iso_val - s0) / (s1 - s0);
        if (tm.is_border(h)) {
            h = tm.opposite(h);
            t = 1 - t;
        }
        auto loc = PMP::locate_on_halfedge(h, t, tm);

        auto vdp = boost::add_vertex(al_graph);
        al_graph[vdp] = loc;
        if (!tm.is_border(h)) {
            auto &&[itm, new_insertion] =
                al_edge_map.insert(std::pair<HalfedgeIndex, AL_vertex_pair>(
                    h, AL_vertex_pair(vdp, AL_graph::null_vertex())));
            if (!new_insertion) {
                assert(itm->second.second == AL_graph::null_vertex());
                itm->second.second = vdp;
                boost::add_edge(itm->second.first, itm->second.second,
                                al_graph);
            }
        }
        h = tm.opposite(h);
        if (!tm.is_border(h)) {
            auto &&[itm, new_insertion] =
                al_edge_map.insert(std::pair<HalfedgeIndex, AL_vertex_pair>(
                    h, AL_vertex_pair(vdp, AL_graph::null_vertex())));
            if (!new_insertion) {
                assert(itm->second.second == AL_graph::null_vertex());
                itm->second.second = vdp;
                boost::add_edge(itm->second.first, itm->second.second,
                                al_graph);
            }
        }
    }

    Polyline_visitor visitor(al_graph);
    CGAL::split_graph_into_polylines(al_graph, visitor);
    return visitor.polylines;
}
} // namespace ic
