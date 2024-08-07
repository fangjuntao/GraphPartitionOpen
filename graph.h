// graph.h
#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include<set>
/* 
class Graph {  //先处理好文件再读图，就不用存这么多乱七八糟的东西了
public:
    std::unordered_map<int, std::vector<double>> nodes;
    std::unordered_map<int, std::unordered_set<int>> edges;  //vector 邻居
    std::vector<double> feat_avg;
  // std::unordered_map<int, int> orig_id_to_graph_id;
   // std::unordered_map<int, int> graph_id_to_orig_id;
    std::unordered_map<int, int> node_id_to_partID;  //node在graph 里面的id  to  partID
    // std::unordered_map<int, int> node_id_to_clusterID; //node在graph 里面的id  to  clusterID 
    int n , m ; 
    void add_node(int orig_node_id, const std::vector<double>& feat);
    void add_edge(int orig_node1, int orig_node2);
    // void normalize();
}; */



//先处理好数据 id 从 0 开始，没有重边
class Graph 
{           //先处理好文件再读图，就不用存这么多乱七八糟的东西了
public:
    std::vector<std::vector<double>> nodes;     //graphID-- features
    std::vector<std::vector<int>> edges;  //graphID-- 邻居
    std::vector<double> feat_sum;
    std::unordered_map<int, int> orig_id_to_graph_id;
    std::unordered_map<int, int> graph_id_to_orig_id;
    std::unordered_map<int, int> node_id_to_partID;  //node在graph 里面的id  to  partID
    std::unordered_map<int, int> node_id_to_clusterID; //node在graph 里面的id  to  clusterID 
    
    
    int n , m ; 
    // void add_node(int orig_node_id, const std::vector<double>& feat);
    // void add_edge(int orig_node1, int orig_node2);

    // void add_node_process(int orig_node_id, const std::vector<double>& feat);
    // void add_edge_process(int orig_node1, int orig_node2);
    // void normalize();
};

Graph read_graph(const std::string& data_input);
Graph read_graphProcressData(const std::string& data_input);
void get_need_file(std::string path, std::vector<std::string>& input_files, std:: string ext);