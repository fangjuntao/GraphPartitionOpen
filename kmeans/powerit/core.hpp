#ifndef POWERIT_CORE_HPP
#define POWERIT_CORE_HPP
 
#include <vector>
#include <cmath>
#include <algorithm>
#include "../aarand/aarand.hpp"
 
namespace powerit {
 
struct Options {
    int iterations = 500;
 
    double tolerance = 0.000001;
 
    int num_threads = 1;
};
 
template<typename Data_>
Data_ normalize(int ndim, Data_* x) {
    Data_ ss = 0;
    for (int d = 0; d < ndim; ++d) {
        ss += x[d] * x[d];
    }
 
    if (ss) {
        ss = std::sqrt(ss);
        for (int d = 0; d < ndim; ++d) {
            x[d] /= ss;
        }
    }
    return ss;
}
template<typename Data_>
struct Result {
    Data_ value;
 
    int iterations;
};
 
template<typename Data_, class Engine_>
void fill_starting_vector(size_t order, Data_* vector, Engine_& engine) {
    while (1) {
        for (size_t d = 1; d < order; d += 2) {
            auto sampled = aarand::standard_normal<Data_>(engine);
            vector[d - 1] = sampled.first;
            vector[d] = sampled.second;
        }
        if (order % 2) {
            vector[order - 1] = aarand::standard_normal<Data_>(engine).first;
        }
        if (normalize(order, vector)) {
            break;
        }
    }
}
 
template<class Multiply_, typename Data_>
Result<Data_> compute_core(size_t order, Multiply_ multiply, Data_* vector, const Options& opt) {
    Result<Data_> stats;
    auto& l2 = stats.value;
    stats.iterations = -1;
    std::vector<Data_> buffer(order);
 
    for (int i = 0; i < opt.iterations; ++i) {
        multiply(buffer, vector);
        l2 = normalize(order, buffer.data());
 
        // Assuming convergence if the vector did not change much from the last iteration.
        Data_ err = 0;
        for (size_t d = 0; d < order; ++d) {
            Data_ diff = buffer[d] - vector[d];
            err += diff * diff;
        }
        if (std::sqrt(err) < opt.tolerance) {
            stats.iterations = i + 1;
            break;
        }
 
        std::copy(buffer.begin(), buffer.end(), vector);
    }
 
    return stats;
} 
 
}
 
#endif
