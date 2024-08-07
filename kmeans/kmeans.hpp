#ifndef KMEANS_KMEANS_HPP
#define KMEANS_KMEANS_HPP
 
#include "Details.hpp"
#include "Refine.hpp"
#include "Initialize.hpp"
#include "MockMatrix.hpp"
 
#include "InitializeKmeanspp.hpp"
#include "InitializeRandom.hpp"
#include "InitializePcaPartition.hpp"
#include "InitializeNone.hpp"
 
#include "RefineHartiganWong.hpp"
#include "RefineLloyd.hpp"
#include "RefineMiniBatch.hpp"
 
#include "compute_wcss.hpp"
 
namespace kmeans {
 
template<class Matrix_, typename Cluster_, typename Float_>
Details<typename Matrix_::index_type> compute(
    const Matrix_& data, 
    const Initialize<Matrix_, Cluster_, Float_>& initialize, 
    const Refine<Matrix_, Cluster_, Float_>& refine,
    Cluster_ num_centers,
    Float_* centers,
    Cluster_* clusters)
{
    auto actual_centers = initialize.run(data, num_centers, centers);
    auto output = refine.run(data, actual_centers, centers, clusters);
    output.sizes.resize(num_centers); // restoring the full size.
    return output;
}
 
template<typename Cluster_, typename Float_, typename Index_>
struct Results {
    template<typename Dim_>
    Results(Dim_ num_dimensions, Index_ num_observations, Cluster_ num_centers) : 
        centers(num_dimensions * num_centers), clusters(num_observations) {}
 
    Results() = default;
    std::vector<Cluster_> clusters;
 
    std::vector<Float_> centers;
 
    Details<Index_> details;
};
 
template<class Matrix_, typename Cluster_, typename Float_>
Results<Cluster_, Float_, typename Matrix_::index_type> compute(
    const Matrix_& data, 
    const Initialize<Matrix_, Cluster_, Float_>& initialize, 
    const Refine<Matrix_, Cluster_, Float_>& refine,
    Cluster_ num_centers)
{
    Results<Cluster_, Float_, typename Matrix_::index_type> output;
    output.clusters.resize(data.num_observations());
    output.centers.resize(static_cast<size_t>(num_centers) * static_cast<size_t>(data.num_dimensions()));
    output.details = compute(data, initialize, refine, num_centers, output.centers.data(), output.clusters.data());
    return output;
}
 
}
 
#endif