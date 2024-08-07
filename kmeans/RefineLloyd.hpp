#ifndef KMEANS_LLOYD_HPP
#define KMEANS_LLOYD_HPP
 
#include <vector>
#include <algorithm>
 
#include "Refine.hpp"
#include "Details.hpp"
#include "QuickSearch.hpp"
#include "is_edge_case.hpp"
#include "compute_centroids.hpp"
#include "parallelize.hpp"
 
namespace kmeans {
 
struct RefineLloydOptions {
    int max_iterations = 10;
 
    int num_threads = 1;
};
 
template<typename Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class RefineLloyd : public Refine<Matrix_, Cluster_, Float_> {
private:
    RefineLloydOptions my_options;
 
    typedef typename Matrix_::index_type Index_;
 
public:
    RefineLloyd(RefineLloydOptions options) : my_options(std::move(options)) {}
 
    RefineLloyd() = default;
 
public:
    RefineLloydOptions& get_options() {
        return my_options;
    }
 
public:
    Details<Index_> run(const Matrix_& data, Cluster_ ncenters, Float_* centers, Cluster_* clusters) const {
        auto nobs = data.num_observations();
        if (internal::is_edge_case(nobs, ncenters)) {
            return internal::process_edge_case(data, ncenters, centers, clusters);
        }
 
        int iter = 0, status = 0;
        std::vector<Index_> sizes(ncenters);
        std::vector<Cluster_> copy(nobs);
        auto ndim = data.num_dimensions();
        internal::QuickSearch<Float_, Cluster_, decltype(ndim)> index;
 
        for (iter = 1; iter <= my_options.max_iterations; ++iter) {
            index.reset(ndim, ncenters, centers);
            internal::parallelize(nobs, my_options.num_threads, [&](int, Index_ start, Index_ length) {
                auto work = data.create_workspace(start, length);
                for (Index_ obs = start, end = start + length; obs < end; ++obs) {
                    auto dptr = data.get_observation(work);
                    copy[obs] = index.find(dptr); 
                }
            });
 
            // Checking if it already converged.
            bool updated = false;
            for (Index_ obs = 0; obs < nobs; ++obs) {
                if (copy[obs] != clusters[obs]) {
                    updated = true;
                    break;
                }
            }
            if (!updated) {
                break;
            }
            std::copy(copy.begin(), copy.end(), clusters);
 
            std::fill(sizes.begin(), sizes.end(), 0);
            for (Index_ obs = 0; obs < nobs; ++obs) {
                ++sizes[clusters[obs]];
            }
            internal::compute_centroids(data, ncenters, centers, clusters, sizes);
        }
 
        if (iter == my_options.max_iterations + 1) {
            status = 2;
        }
 
        return Details<Index_>(std::move(sizes), iter, status);
    }
};
 
}
 
#endif
