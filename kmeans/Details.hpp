#ifndef KMEANS_DETAILS_HPP
#define KMEANS_DETAILS_HPP
 
#include <vector>
 
namespace kmeans {
 
template<typename Index_ = int>
struct Details {
    Details() = default;
 
    Details(int iterations, int status) : sizes(0), iterations(iterations), status(status) {}
 
    Details(std::vector<Index_> sizes, int iterations, int status) : sizes(std::move(sizes)), iterations(iterations), status(status) {} 
    std::vector<Index_> sizes;
 
    int iterations = 0;
 
    int status = 0;
};
 
}
 
#endif
 