#ifndef KMEANS_REFINE_HPP
#define KMEANS_REFINE_HPP
 
#include "Details.hpp"
#include "SimpleMatrix.hpp"
 
namespace kmeans {
 
template<typename Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class Refine {
public:
    virtual ~Refine() = default;
 
    virtual Details<typename Matrix_::index_type> run(const Matrix_& data, Cluster_ num_centers, Float_* centers, Cluster_* clusters) const = 0;
};
 
}
 
#endif