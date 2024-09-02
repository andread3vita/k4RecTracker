#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <span>
#include <utility>
#include <vector>
#include <torch/torch.h>
#include <nanoflann/nanoflann.hpp>

/**
 * @struct point3w
 * A structure representing a 3D point with an associated weight.
 *
 * This structure is used to store the x, y, z coordinates of a point
 * along with an additional float value which could represent a weight
 * or beta value.
 */
struct point3w {
    double x, y, z, beta;
};

/**
 * @dbscan Convenience function to perform DBSCAN on point3w data.
 *
 * This function wraps the general dbscan function for use with
 * point3w data. It initializes the necessary adaptor and performs
 * the clustering.
 *
 * @param data A span of point3w objects representing the 3D points to be clustered.
 * @param epsilon The distance threshold for clustering (epsilon).
 * @param min_pts The minimum number of points required to form a cluster.
 * @return A vector of clusters, where each cluster is a vector of point indices.
 */
auto dbscan(const std::span<const point3w>& data, double epsilon, int min_pts) -> std::vector<std::vector<size_t>>;