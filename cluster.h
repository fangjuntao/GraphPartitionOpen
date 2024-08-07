// cluster.h
#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
class Cluster {   // 存储partition 信息
public:
    int count;  //节点的数量
    std::vector<double> vector;  //特征
    std::unordered_map<int, int> point;  // 一个点是否在该分区内   为0  代表不在， 大于 0 代表在 （1）
    std::unordered_set<int> nodeSet;    //点集
    long long  innerEdgeNum;  //内部边数
    int outerEdgeNum;  //外部边数
    long long  edgeSum;      //总边数

    Cluster(int count = 0, const std::vector<double>& vector = std::vector<double>(), long long innerEdgeNum = 0, int outerEdgeNum = 0, long long edgeSum = 0);


};
