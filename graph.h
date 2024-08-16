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




class Graph 
{           
public:
    std::vector<std::vector<double>> nodes;     //graphID-- features
    std::vector<std::vector<int>> edges;  //graphID
    std::vector<double> feat_sum;
    std::unordered_map<int, int> orig_id_to_graph_id;
    std::unordered_map<int, int> graph_id_to_orig_id;
    std::unordered_map<int, int> node_id_to_partID;  
    std::unordered_map<int, int> node_id_to_clusterID; 
    
    
    int n , m ; 

};

Graph read_graph(const std::string& data_input);
Graph read_graphProcressData(const std::string& data_input);
void get_need_file(std::string path, std::vector<std::string>& input_files, std:: string ext);
