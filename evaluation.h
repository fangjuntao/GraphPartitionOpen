// evaluation.h
#pragma once
#include <vector>
#include "cluster.h"
#include "graph.h"
using FiveTuple = std::tuple<bool, int, int, int, int>;




std::vector<double> computeObject(std::vector<Cluster>& clusters, std::vector<int>& nodeNumList, std:: vector<std::vector<double>>& featuresList, std::vector<int>& innerEdgeList ,  std::vector<int>& outerEdgeList, const Graph& g, double count_optimal, const std::vector<double>& vector_optimal,  const int & k );


std::tuple<int, std::vector<std::unordered_map<std::string, std::vector<std::pair<int, int>>>>, std::unordered_map<int, std::string>>
calculateGain(const Graph& g, const std::vector<std::vector<int>>& cluster_node_list, std::unordered_map<int, std::vector<int>>& node2neighborsPart, int npart);



void getIndex(const std::string& key, int& index_1, int& index_2);


FiveTuple IfSatisfiedContrainst(const Graph & g , const int& node, const int & index_1, const int & index_2, const std::unordered_map<int, std::vector<int>> & node2neighborsPart, std::vector<int> nodeNumListTmp, std::vector<std::vector<double>> featuresListTmp, std::vector<int> innerEdgeListTmp, std::vector<int> outerEdgeListTmp, double imbalance1, double imbalance2, double imbalance3, double imbalance4, int count_optimal, const std::vector<double>& vector_optimal) ;//O(1)


void saveResult(std::vector<Cluster>& result, const Graph& g, const int& feacuter_size, const int& k) ;



std::tuple<int, std::unordered_map<int, std::pair<int,int>>,std::unordered_map<int, int>>
calculateGain(const Graph& g, std::unordered_map<int, std::vector<int>>& node2neighborsPart, const int& npart);



void saveResultForIndicator(std::vector<Cluster>& result, Graph& g, const int& feacuter_size, const int& k) ;


