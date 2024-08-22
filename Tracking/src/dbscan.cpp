#include "dbscan.hpp"

inline auto get_pt(const point3w& p, std::size_t dim) {
    if (dim == 0) return p.x;
    if (dim == 1) return p.y;
    return p.z;
}

template<typename Point>
struct adaptor {
    const std::span<const Point>& points;
    adaptor(const std::span<const Point>& pts) : points(pts) { } // Cambiato da 'points' a 'pts'

    inline std::size_t kdtree_get_point_count() const { return points.size(); }

    inline float kdtree_get_pt(const std::size_t idx, const std::size_t dim) const {
        return get_pt(points[idx], dim);
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

    auto const * elem_ptr(const std::size_t idx) const {
        return &points[idx].x;
    }
};

auto sort_clusters(std::vector<std::vector<size_t>>& clusters) {
    for (auto& cluster : clusters) {
        std::sort(cluster.begin(), cluster.end());
    }
}

template<typename Adaptor>
auto dbscan(const Adaptor& adapt, float epsilon, int min_pts) {
    epsilon *= epsilon;
    using namespace nanoflann;
    using my_kd_tree_t = KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, decltype(adapt)>, decltype(adapt), 3>;

    auto index = my_kd_tree_t(3, adapt, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    const auto n_points = adapt.kdtree_get_point_count();
    auto visited = std::vector<bool>(n_points);
    auto clusters = std::vector<std::vector<size_t>>();
    auto matches = std::vector<std::pair<size_t, float>>();
    auto sub_matches = std::vector<std::pair<size_t, float>>();

    for (size_t i = 0; i < n_points; i++) {
        if (visited[i]) continue;

        index.radiusSearch(adapt.elem_ptr(i), epsilon, matches, SearchParams(32, 0.f, false));
        if (matches.size() < static_cast<size_t>(min_pts)) continue;
        visited[i] = true;

        auto cluster = std::vector({i});

        while (!matches.empty()) {
            auto nb_idx = matches.back().first;
            matches.pop_back();
            if (visited[nb_idx]) continue;
            visited[nb_idx] = true;

            index.radiusSearch(adapt.elem_ptr(nb_idx), epsilon, sub_matches, SearchParams(32, 0.f, false));

            if (sub_matches.size() >= static_cast<size_t>(min_pts)) {
                std::copy(sub_matches.begin(), sub_matches.end(), std::back_inserter(matches));
            }
            cluster.push_back(nb_idx);
        }
        clusters.emplace_back(std::move(cluster));
    }
    sort_clusters(clusters);
    return clusters;
}

auto dbscan(const std::span<const point3w>& data, float epsilon, int min_pts) -> std::vector<std::vector<size_t>> {
    const auto adapt = adaptor<point3w>(data);
    return dbscan(adapt, epsilon, min_pts);
}
