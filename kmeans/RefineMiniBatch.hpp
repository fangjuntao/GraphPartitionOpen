#ifndef KMEANS_REFINE_MINIBATCH_HPP
#define KMEANS_REFINE_MINIBATCH_HPP
 
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <random>
#include <type_traits>
 
#include "aarand/aarand.hpp"
 
#include "Refine.hpp"
#include "Details.hpp"
#include "QuickSearch.hpp"
#include "is_edge_case.hpp"
#include "parallelize.hpp"
 
namespace kmeans {
 
struct RefineMiniBatchOptions {
    int max_iterations = 100;
 
    int batch_size = 500;
 
    double max_change_proportion = 0.01;
 
    int convergence_history = 10;
 
    uint64_t seed = 1234567890u;
 
    int num_threads = 1;
};
 
template<typename Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class RefineMiniBatch : public Refine<Matrix_, Cluster_, Float_> {
public:
    RefineMiniBatch(RefineMiniBatchOptions options) : my_options(std::move(options)) {}
 
    RefineMiniBatch() = default;
 
public:
    RefineMiniBatchOptions& get_options() {
        return my_options;
    }
 
private:
    RefineMiniBatchOptions my_options;
 
public:
    Details<typename Matrix_::index_type> run(const Matrix_& data, Cluster_ ncenters, Float_* centers, Cluster_* clusters) const {
        auto nobs = data.num_observations();
        if (internal::is_edge_case(nobs, ncenters)) {
            return internal::process_edge_case(data, ncenters, centers, clusters);
        }
 
        int iter = 0, status = 0;
        std::vector<uint64_t> total_sampled(ncenters); // holds the number of sampled observations across iterations, so we need a large integer.
        std::vector<Cluster_> previous(nobs);
        typedef decltype(nobs) Index_;
        std::vector<uint64_t> last_changed(ncenters), last_sampled(ncenters); // holds the number of sampled/changed observation for the last few iterations.
 
        Index_ actual_batch_size = nobs;
        typedef typename std::conditional<std::is_signed<Index_>::value, int, unsigned int>::type SafeCompInt; // waiting for C++20's comparison functions...
        if (static_cast<SafeCompInt>(actual_batch_size) > my_options.batch_size) {
            actual_batch_size = my_options.batch_size;
        }
        std::vector<Index_> chosen(actual_batch_size);
        std::mt19937_64 eng(my_options.seed);
 
        auto ndim = data.num_dimensions();
        size_t long_ndim = ndim;
        internal::QuickSearch<Float_, Cluster_, decltype(ndim)> index;
 
        for (iter = 1; iter <= my_options.max_iterations; ++iter) {
            aarand::sample(nobs, actual_batch_size, chosen.data(), eng);
            if (iter > 1) {
                for (auto o : chosen) {
                    previous[o] = clusters[o];
                }
            }
 
            index.reset(ndim, ncenters, centers);
            internal::parallelize(actual_batch_size, my_options.num_threads, [&](int, Index_ start, Index_ length) {
                auto work = data.create_workspace(chosen.data() + start, length);
                for (Index_ s = start, end = start + length; s < end; ++s) {
                    auto ptr = data.get_observation(work);
                    clusters[chosen[s]] = index.find(ptr);
                }
            });
 
            // Updating the means for each cluster.
            auto work = data.create_workspace(chosen.data(), actual_batch_size);
            for (auto o : chosen) {
                const auto c = clusters[o];
                auto& n = total_sampled[c];
                ++n;
 
                Float_ mult = static_cast<Float_>(1)/static_cast<Float_>(n);
                auto ccopy = centers + static_cast<size_t>(c) * long_ndim;
                auto ocopy = data.get_observation(work);
 
                for (decltype(ndim) d = 0; d < ndim; ++d, ++ocopy, ++ccopy) {
                    (*ccopy) += (static_cast<Float_>(*ocopy) - *ccopy) * mult; // cast to ensure consistent precision regardless of Matrix_::data_type.
                }
            }
 
            // Checking for updates.
            if (iter != 1) {
                for (auto o : chosen) {
                    auto p = previous[o];
                    ++(last_sampled[p]);
                    auto c = clusters[o];
                    if (p != c) {
                        ++(last_sampled[c]);
                        ++(last_changed[p]);
                        ++(last_changed[c]);
                    }
                }
 
                if (iter % my_options.convergence_history == 1) {
                    bool too_many_changes = false;
                    for (Cluster_ c = 0; c < ncenters; ++c) {
                        if (static_cast<double>(last_changed[c]) >= static_cast<double>(last_sampled[c]) * my_options.max_change_proportion) {
                            too_many_changes = true;
                            break;
                        }
                    }
 
                    if (!too_many_changes) {
                        break;
                    }
                    std::fill(last_sampled.begin(), last_sampled.end(), 0);
                    std::fill(last_changed.begin(), last_changed.end(), 0);
                }
            }
        }
 
        if (iter == my_options.max_iterations + 1) {
            status = 2;
        }
 
        // Run through all observations to make sure they have the latest cluster assignments.
        index.reset(ndim, ncenters, centers);
        internal::parallelize(nobs, my_options.num_threads, [&](int, Index_ start, Index_ length) {
            auto work = data.create_workspace(start, length);
            for (Index_ s = start, end = start + length; s < end; ++s) {
                auto ptr = data.get_observation(work);
                clusters[s] = index.find(ptr);
            }
        });
 
        std::vector<Index_> cluster_sizes(ncenters);
        for (Index_ o = 0; o < nobs; ++o) {
            ++cluster_sizes[clusters[o]];
        }
 
        internal::compute_centroids(data, ncenters, centers, clusters, cluster_sizes);
        return Details<Index_>(std::move(cluster_sizes), iter, status);
    }
};
 
}
 
#endif
