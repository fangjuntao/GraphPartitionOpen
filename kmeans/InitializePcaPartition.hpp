#ifndef KMEANS_INITIALIZE_PCA_PARTITION_HPP
#define KMEANS_INITIALIZE_PCA_PARTITION_HPP
 
#include <vector>
#include <algorithm>
#include <numeric>
#include <queue>
#include <random>
#include <cstdint>
 
#include "aarand/aarand.hpp"
#include "powerit/powerit.hpp"
 
#include "Initialize.hpp"
#include "compute_centroids.hpp"
 
namespace kmeans {
 
struct InitializePcaPartitionOptions {
    powerit::Options power_iteration_options;
 
    double size_adjustment = 1;
 
    uint64_t seed = 6523u;
};
 
namespace InitializePcaPartition_internal {
 
template<typename Float_>
struct Workspace {
    Workspace(size_t ndim) : pc(ndim), delta(ndim), cov(ndim * ndim) {}
    std::vector<Float_> pc;
    std::vector<Float_> delta;
    std::vector<Float_> cov;
};
 
template<class Matrix_, typename Float_, class Engine_>
void compute_pc1(
    const Matrix_& data, 
    const std::vector<typename Matrix_::index_type>& chosen, 
    const Float_* center, 
    Engine_& eng, 
    Workspace<Float_>& work,
    const powerit::Options& power_opts)
{
    auto ndim = data.num_dimensions();
    size_t long_ndim = ndim;
 
    auto matwork = data.create_workspace(chosen.data(), chosen.size());
    for (size_t i = 0, end = chosen.size(); i < end; ++i) {
        auto dptr = data.get_observation(matwork);
        for (decltype(ndim) j = 0; j < ndim; ++j) {
            work.delta[j] = static_cast<Float_>(dptr[j]) - center[j]; // cast to ensure consistent precision regardless of Data_.
        }
 
        size_t offset = 0;
        for (decltype(ndim) j = 0; j < ndim; ++j, offset += long_ndim) {
            size_t copy_offset = offset;
            for (decltype(ndim) k = 0; k <= j; ++k, ++copy_offset) {
                work.cov[copy_offset] += work.delta[j] * work.delta[k];
            }
        }
    }
 
    // Filling in the other side of the matrix, to enable cache-efficient multiplication.
    size_t src_offset = 0;
    for (decltype(ndim) j = 0; j < ndim; ++j, src_offset += long_ndim) {
        size_t dest_offset = j;
        size_t src_offset_copy = src_offset;
        for (decltype(ndim) k = 0; k <= j; ++k, ++src_offset_copy, dest_offset += long_ndim) {
            work.cov[dest_offset] = work.cov[src_offset_copy];
        }
    }
 
    powerit::compute(ndim, work.cov.data(), /* row_major = */ true, work.pc.data(), eng, power_opts);
} 
 
template<class Matrix_, typename Float_>
void compute_center(const Matrix_& data, const std::vector<typename Matrix_::index_type>& chosen, Float_* center) {
    auto ndim = data.num_dimensions();
    std::fill_n(center, ndim, 0);
 
    auto work = data.create_workspace(chosen.data(), chosen.size());
    for (size_t i = 0, end = chosen.size(); i < end; ++i) {
        auto dptr = data.get_observation(work);
        for (decltype(ndim) d = 0; d < ndim; ++d) {
            center[d] += static_cast<Float_>(dptr[d]); // cast to ensure consistent precision regardless of Data_.
        }
    }
 
    for (decltype(ndim) d = 0; d < ndim; ++d) {
        center[d] /= chosen.size();
    }
}
 
template<class Matrix_, typename Float_>
Float_ update_center_and_mrse(const Matrix_& data, const std::vector<typename Matrix_::index_type>& chosen, Float_* center) {
    compute_center(data, chosen, center);
 
    auto ndim = data.num_dimensions();
    Float_ mrse = 0;
    auto work = data.create_workspace(chosen.data(), chosen.size());
    for (size_t i = 0, end = chosen.size(); i < end; ++i) {
        auto dptr = data.get_observation(work);
        for (decltype(ndim) d = 0; d < ndim; ++d) {
            Float_ delta = static_cast<Float_>(dptr[d]) - center[d]; // cast to ensure consistent precision regardless of Data_.
            mrse += delta * delta;
        }
    }
 
    return mrse / chosen.size();
}
 
}
template<typename Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class InitializePcaPartition : public Initialize<Matrix_, Cluster_, Float_> {
public:
    InitializePcaPartition(InitializePcaPartitionOptions options) : my_options(std::move(options)) {}
 
    InitializePcaPartition() = default;
 
private:
    InitializePcaPartitionOptions my_options;
 
public:
    InitializePcaPartitionOptions& get_options() {
        return my_options;
    }
 
public:
    Cluster_ run(const Matrix_& data, Cluster_ ncenters, Float_* centers) const {
        auto nobs = data.num_observations();
        if (nobs == 0) {
            return 0;
        }
 
        std::mt19937_64 rng(my_options.seed);
        std::priority_queue<std::pair<Float_, Cluster_> > mrse;
        std::vector<std::vector<typename Matrix_::index_type> > assignments(ncenters);
 
        auto ndim = data.num_dimensions();
        InitializePcaPartition_internal::Workspace<Float_> power_work(ndim);
 
        // Setting up the zero'th cluster. (No need to actually compute the
        // MRSE at this point, as there's nothing to compare it to.)
        internal::compute_centroid(data, centers);
        assignments[0].resize(nobs);
        std::iota(assignments.front().begin(), assignments.front().end(), 0);
        std::vector<typename Matrix_::index_type> replace_assignments;
 
        for (Cluster_ cluster = 1; cluster < ncenters; ++cluster) {
            Cluster_ worst_cluster = 0;
            if (mrse.size()) {
                worst_cluster = mrse.top().second;
                mrse.pop();
            }
 
            // Extracting the principal component for this bad boy.
            auto worst_center = centers + static_cast<size_t>(worst_cluster) * static_cast<size_t>(ndim); // cast to avoid overflow.
            auto& worst_assignments = assignments[worst_cluster];
            InitializePcaPartition_internal::compute_pc1(data, worst_assignments, worst_center, rng, power_work, my_options.power_iteration_options);
            const auto& pc1 = power_work.pc;
 
            // Projecting all points in this cluster along PC1. The center lies
            // at zero, so everything positive (on one side of the hyperplane
            // orthogonal to PC1 and passing through the center) gets bumped to
            // the next cluster.
            auto& new_assignments = assignments[cluster];
            replace_assignments.clear();
 
            size_t num_in_cluster = worst_assignments.size();
            auto work = data.create_workspace(worst_assignments.data(), num_in_cluster);
            for (auto i : worst_assignments) {
                auto dptr = data.get_observation(work);
 
                Float_ proj = 0;
                for (decltype(ndim) d = 0; d < ndim; ++d, ++dptr) {
                    proj += (static_cast<Float_>(*dptr) - worst_center[d]) * pc1[d]; // cast for consistent precision regardless of Data_.
                }
 
                if (proj > 0) {
                    new_assignments.push_back(i);
                } else {
                    replace_assignments.push_back(i);
                }
            }
 
            // If one or the other is empty, then this entire procedure short
            // circuits as all future iterations will just re-select this
            // cluster (which won't get partitioned properly anyway). In the
            // bigger picture, the quick exit out of the iterations is correct
            // as we should only fail to partition in this manner if all points
            // within each remaining cluster are identical.
            if (new_assignments.empty() || replace_assignments.empty()) {
                return cluster;
            }
 
            // Computing centers and MRSE.
            auto new_center = centers + static_cast<size_t>(cluster) * static_cast<size_t>(ndim); // cast to avoid overflow.
            auto new_mrse = InitializePcaPartition_internal::update_center_and_mrse(data, new_assignments, new_center);
            mrse.emplace(new_mrse, cluster);
 
            auto replace_mrse = InitializePcaPartition_internal::update_center_and_mrse(data, replace_assignments, worst_center);
            mrse.emplace(replace_mrse, worst_cluster);
 
            worst_assignments.swap(replace_assignments);
        }
 
        return ncenters;
    }
};
 
}
 
#endif