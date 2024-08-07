#ifndef KMEANS_INITIALIZE_HPP
#define KMEANS_INITIALIZE_HPP
 
#include "SimpleMatrix.hpp"
 
namespace kmeans {
 
template<class Matrix_ = SimpleMatrix<double, int>, typename Cluster_ = int, typename Float_ = double>
class Initialize {
public:
    virtual ~Initialize() = default;
    virtual Cluster_ run(const Matrix_& data, Cluster_ num_centers, Float_* centers) const = 0;
};
 
}
 
#endif