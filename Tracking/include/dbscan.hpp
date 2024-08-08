#pragma once

#include <cassert>
#include <cstddef>
#include <span>
#include <vector>
#include <cstdlib>

struct point3w {
    float x, y, z, weight;
};

auto dbscan(const std::span<const point3w>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>;
