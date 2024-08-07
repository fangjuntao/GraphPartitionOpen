#ifndef KMEANS_COMPUTE_WCSS_HPP
#define KMEANS_COMPUTE_WCSS_HPP
 
#include <algorithm>
#include "SimpleMatrix.hpp"
 
namespace kmeans {
 
template<class Matrix_, typename Cluster_, typename Float_>
void compute_wcss(const Matrix_& data, Cluster_ ncenters, const Float_* centers, const Cluster_* clusters, Float_* wcss) {
    auto nobs = data.num_observations();
    auto ndim = data.num_dimensions();
    size_t long_ndim = ndim;
    std::fill_n(wcss, ncenters, 0);
 
    auto work = data.create_workspace(static_cast<decltype(nobs)>(0), nobs);
    for (decltype(nobs) obs = 0; obs < nobs; ++obs) {
        auto cen = clusters[obs];
        auto curcenter = centers + static_cast<size_t>(cen) * long_ndim;
        auto& curwcss = wcss[cen];
 
        auto curdata = data.get_observation(work);
        for (decltype(ndim) dim = 0; dim < ndim; ++dim, ++curcenter, ++curdata) {
            Float_ delta = static_cast<Float_>(*curdata) - *curcenter; // cast for consistent precision regardless of Data_.
            curwcss += delta * delta;
        }
    }
}
 
}
 
#endif