// cluster.cpp
#include "cluster.h"

Cluster::Cluster(int count, const std::vector<double>& vector, long long innerEdgeNum, int outerEdgeNum, long long edgeSum)
    : count(count), vector(vector), innerEdgeNum(innerEdgeNum), outerEdgeNum(outerEdgeNum), edgeSum(edgeSum) {
}