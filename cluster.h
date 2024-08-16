// cluster.h
#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
class Cluster {  
public:
    int count;  
    std::vector<double> vector;  
    std::unordered_map<int, int> point;  
    std::unordered_set<int> nodeSet;    
    long long  innerEdgeNum;  
    int outerEdgeNum;  
    long long  edgeSum;     

    Cluster(int count = 0, const std::vector<double>& vector = std::vector<double>(), long long innerEdgeNum = 0, int outerEdgeNum = 0, long long edgeSum = 0);


};
