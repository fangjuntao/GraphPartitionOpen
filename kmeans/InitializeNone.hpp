#ifndef KMEANS_INITIALIZE_NONE_HPP
#define KMEANS_INITIALIZE_NONE_HPP 
 
#include "Initialize.hpp"
#include <algorithm>
 
namespace kmeans {
 
template<class Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class InitializeNone : public Initialize<Matrix_, Cluster_, Float_> { 
public:
    Cluster_ run(const Matrix_& matrix, Cluster_ ncenters, Float_*) const {
        return std::min(matrix.num_observations(), static_cast<typename Matrix_::index_type>(ncenters));
    }
};
 
}
 
#endif