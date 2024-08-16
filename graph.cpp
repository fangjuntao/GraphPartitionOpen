// graph.cpp
#include "graph.h"
#include <dirent.h>
#include <algorithm>


std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}




void get_need_file(std::string path, std::vector<std::string>& input_files, std:: string ext)
{

        DIR* dir = opendir(path.c_str());
        if (dir) 
        {
            struct dirent* entry;
            while ((entry = readdir(dir)) != nullptr) 
            {
                std::string file_name = entry->d_name;
                if (file_name.find(ext) != std::string::npos) 
                {
                    input_files.push_back(path+"/"+file_name);
                }        
            }
            closedir(dir);
        } 
        else 
        {
            std::cout << "Error: Unable to open directory" << std::endl;
            
        }

  

}



Graph read_graphProcressData(const std::string& data_input) 
{
    Graph g;
    g.m = 0 ;
    g.n = 0 ; 
    std::unordered_map<int, int>&  name2id  =  g.orig_id_to_graph_id;
    std::unordered_map<int, int>& id2name  = g.graph_id_to_orig_id ;   
    std::unordered_map<int, std::unordered_set<int>> graph;   
    std::vector<std::string> paths = split(data_input, ',');

    std::string line;

    std::vector<std::string> input_files;
    

    // Read edges

    std::string edgesPath = paths[0]  ;

    std::string need_extension = "csv";
    get_need_file(edgesPath, input_files, need_extension);

    std::sort(input_files.begin(), input_files.end());
    // std::set<int> nodeSet ; 

    int x , y ;
    for (const auto& file : input_files) 
    {

        std::ifstream edge_file(file);
        
        while (std::getline(edge_file, line)) 
        {
            std::vector<std::string> tokens = split(line, ',');
            int orig_node1 = std::stoi(tokens[0]);
            int orig_node2 = std::stoi(tokens[1]);
            if (name2id.find(orig_node1) == name2id.end())
                name2id[orig_node1] = g.n, id2name[g.n] = orig_node1, g.n++;
            if (name2id.find(orig_node2) == name2id.end())
                name2id[orig_node2] = g.n, id2name[g.n] = orig_node2, g.n++;
            
            
            x = name2id[orig_node1];
            y = name2id[orig_node2];
            graph[x].insert(y);
            graph[y].insert(x);


        }
    }
    g.edges.resize(g.n);
    for (int x = 0; x < g.n; x++)
        for (int y : graph[x])
            g.edges[x].push_back(y), g.m++;
    
    
    
    graph.clear();

    g.m /= 2 ;




    std::cout<<"read edges over!"<<std::endl;



    // Read attribute features
    g.nodes.resize(g.n);
    std::string FeaturesPath = paths[1]  ;
    input_files.clear();
    get_need_file(FeaturesPath, input_files, need_extension);

    std::cout<<"input_files:"<<input_files[0]<<std::endl;
    std::sort(input_files.begin(), input_files.end());

    for (const auto& file : input_files) 
    {

        std::ifstream attr_file(file);
    
        while (std::getline(attr_file, line)) 
        {
            std::vector<std::string> tokens = split(line, ',');
            int orig_node_id = std::stoi(tokens[0]);
            std::vector<double> feat(tokens.size() - 1);
            for (size_t i = 1; i < tokens.size(); ++i) 
            {
                feat[i - 1] = std::stof(tokens[i]);
            }
            g.nodes[ g.orig_id_to_graph_id[orig_node_id] ] = feat;

        }


    }

    std::cout<<"read features over!"<<std::endl;






    if(g.nodes.size()<= 0)
    {
        std::cout<<"graph.cpp line 180 : read graph wrong !"<<std::endl;
    }
    else
    { 

        g.feat_sum.resize(g.nodes[0].size(),0.0);
        for (int i = 0 ; i < g.n ; i++)  
        {
            
            for( int j = 0 ; j < g.nodes[i].size(); j++ )
            {
                g.feat_sum[j] += g.nodes[i][j];
            }
                
        }


    }


    return g;
}





Graph read_graph(const std::string &data_input)
{




    Graph g;
    g.m = 0 ;
    g.n = 0 ; 
    std::unordered_map<int, int>&  name2id  =  g.orig_id_to_graph_id;
    std::unordered_map<int, int>& id2name  = g.graph_id_to_orig_id ;   
    std::unordered_map<int, std::unordered_set<int>> graph;   //辅助读图
    std::vector<std::string> paths = split(data_input, ',');

    std::string line;

    std::vector<std::string> input_files;
    

    // Read edges

    std::string edgesPath = paths[0]  ;

    std::string need_extension = "csv";
    get_need_file(edgesPath, input_files, need_extension);

    std::sort(input_files.begin(), input_files.end());
    // std::set<int> nodeSet ; 

    int x , y ;
    for (const auto& file : input_files) 
    {

        std::ifstream edge_file(file);
        
        while (std::getline(edge_file, line)) 
        {
            std::vector<std::string> tokens = split(line, ',');
            int orig_node1 = std::stoi(tokens[0]);
            int orig_node2 = std::stoi(tokens[1]);
            if (name2id.find(orig_node1) == name2id.end())
                name2id[orig_node1] = g.n, id2name[g.n] = orig_node1, g.n++;
            if (name2id.find(orig_node2) == name2id.end())
                name2id[orig_node2] = g.n, id2name[g.n] = orig_node2, g.n++;
            
            
            x = name2id[orig_node1];
            y = name2id[orig_node2];
            graph[x].insert(y);
            graph[y].insert(x);


        }
    }
    g.edges.resize(g.n);
    for (int x = 0; x < g.n; x++)
        for (int y : graph[x])
            g.edges[x].push_back(y), g.m++;
    
    
    
    graph.clear();

    g.m /= 2 ;

    // g.n   =   nodeSet.size(); 


    std::cout<<"read edges over!"<<std::endl;



    // Read attribute features
    g.nodes.resize(g.n);
    std::string FeaturesPath = paths[1]  ;
    input_files.clear();
    get_need_file(FeaturesPath, input_files, need_extension);

    std::cout<<"input_files:"<<input_files[0]<<std::endl;
    std::sort(input_files.begin(), input_files.end());

    for (const auto& file : input_files) 
    {

        std::ifstream attr_file(file);
    
        while (std::getline(attr_file, line)) 
        {
            std::vector<std::string> tokens = split(line, ',');
            int orig_node_id = std::stoi(tokens[0]);
            std::vector<double> feat(tokens.size() - 1);
            for (size_t i = 1; i < tokens.size(); ++i) 
            {
                feat[i - 1] = std::stof(tokens[i]);
            }
            g.nodes[ g.orig_id_to_graph_id[orig_node_id] ] = feat;
            // g.add_node(orig_node_id, feat);
            // g.n++;
        }


    }

    std::cout<<"read features over!"<<std::endl;






    if(g.nodes.size()<= 0)
    {
        std::cout<<"graph.cpp line 180 : read graph wrong !"<<std::endl;
    }
    else
    { 

        g.feat_sum.resize(g.nodes[0].size(),0.0);
        for (int i = 0 ; i < g.n ; i++)  
        {
            
            for( int j = 0 ; j < g.nodes[i].size(); j++ )
            {
                g.feat_sum[j] += g.nodes[i][j];
            }
                
        }


    }


    return g;






}
