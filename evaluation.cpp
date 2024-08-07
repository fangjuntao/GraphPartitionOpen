// evaluation.cpp
#include "evaluation.h"
#include <algorithm>
#include <cmath>



//O(m)  graph  edge   int -> 点数可以到20亿
std::vector<double> computeObject(std::vector<Cluster>& clusters, std::vector<int>& nodeNumList, std:: vector<std::vector<double>>& featuresList, std::vector<int>& innerEdgeList ,  std::vector<int>& outerEdgeList, const Graph& g, double count_optimal, const std::vector<double>& vector_optimal,  const int & k )
{
    

    int  min_out = 0;  // partition 中的最少外边数
    int min_count = 0;
    long long min_innerEdge = 0;
    long long max_innerEdge = 0  ;
    int max_out = 0 ;
    int max_count =  0 ;
    long  long  out_   =  0 ;
    long  long  Edge_ =  0 ;
    nodeNumList.clear();
    featuresList.clear();
    innerEdgeList.clear();  
    outerEdgeList.clear();
    for (int partID = 0; partID < k; ++partID) 
    {
        
 
        if (clusters[partID].nodeSet.size() != static_cast<size_t>(clusters[partID].count)) 
        {
            std::cout << "-------------------------------line 256 wrong----------------------------------" << std::endl;
        }

        for (int point : clusters[partID].nodeSet) 
        {

            const std::vector<int>& neighbors = g.edges[point];
            int innerNum = 0;
            for (int neighbor : neighbors) 
            {
                if (clusters[partID].nodeSet.count(neighbor) > 0)   //O（1） hash set
                {
                    ++innerNum;
                }
            }
            clusters[partID].innerEdgeNum += innerNum;

        }

        clusters[partID].outerEdgeNum = clusters[partID].edgeSum - clusters[partID].innerEdgeNum;
        clusters[partID].innerEdgeNum /= 2;

        if (partID == 0) 
        {
            min_out = clusters[partID].outerEdgeNum;
            min_count = clusters[partID].count;
            min_innerEdge = clusters[partID].innerEdgeNum;
        } 
        else 
        {
            min_out = std::min(min_out, clusters[partID].outerEdgeNum);
            min_count = std::min(min_count, clusters[partID].count);
            min_innerEdge = std::min(min_innerEdge, clusters[partID].innerEdgeNum);
        }

        
        max_out  =  std::max(max_out,  clusters[partID].outerEdgeNum );
        max_count =  std::max(max_count,  clusters[partID].count) ;

        max_innerEdge = std::max(max_innerEdge, clusters[partID].innerEdgeNum );
        out_  +=   clusters[partID].outerEdgeNum;  
        Edge_ +=  clusters[partID].edgeSum;       
        nodeNumList.push_back(clusters[partID].count);
        innerEdgeList.push_back( static_cast<int>(clusters[partID].innerEdgeNum));     // 这里需要确保 clusters[partID].innerEdgeNum 的值小于20 亿，类型转换才不会出错 
        outerEdgeList.push_back(clusters[partID].outerEdgeNum);
        std::vector<double> scaled_vector = clusters[partID].vector;
        int count = clusters[partID].count;
        std::transform(scaled_vector.begin(), scaled_vector.end(), scaled_vector.begin(), [count](double val) { return val * count; });
        featuresList.push_back(scaled_vector);

        
    }
    
    
    out_  /=  2;
    Edge_ /=  2;
    // double innerdge_optimal = static_cast<double>(Edge_) / k;
    
    std::vector<double> ans;

    double max_abs_diff_divided_1 = 0.0;
    std::vector<double> max_abs_diff_divided_2(vector_optimal.size(), 0.0);
    double max_abs_diff_divided_3 = 0.0;
    double min_max_ratio2 = 0.0;




    for (const Cluster& cluster : clusters) 
    {
        int count = cluster.count;
        const std::vector<double>& vector = cluster.vector;
        // int innerEdgeNum =  static_cast<int>(cluster.innerEdgeNum);
        // int outerEdgeNum =  cluster.outerEdgeNum;

        max_abs_diff_divided_1 = std::max(max_abs_diff_divided_1, std::abs(static_cast<double>(count) - count_optimal) / count_optimal);

        for (size_t i = 0; i < vector.size(); ++i)
        {
            max_abs_diff_divided_2[i] = std::max(max_abs_diff_divided_2[i], std::abs(vector[i] * count - vector_optimal[i]) / vector_optimal[i]);
        }

        //max_abs_diff_divided_3 = std::max(max_abs_diff_divided_3, std::abs(static_cast<double>(innerEdgeNum) - innerdge_optimal) / innerdge_optimal);

    }
    
    
    
    max_abs_diff_divided_3  =  (1.0* min_innerEdge) / max_innerEdge ;
    min_max_ratio2  =  ( 1.0*min_out ) / max_out ;
    ans.push_back( max_abs_diff_divided_1 );
    ans.insert( ans.end(), max_abs_diff_divided_2.begin(), max_abs_diff_divided_2.end());
    ans.push_back( max_abs_diff_divided_3 );
    ans.push_back( min_max_ratio2 );
    ans.push_back((1.0* out_ )/ Edge_);

    std::cout<<"溢出效应："<< (1.0* out_ )/ Edge_ << std::endl;
    std::cout<<"溢出边的数量："<< out_<<std::endl;
    std::cout<<"总边的数量："<< Edge_<<std::endl;   
    std::cout<<"纯度："<< 1- (1.0* out_ )/ Edge_ << std::endl;   
    std::cout<<"溢出边的比例："<< ( 1.0*min_out ) / max_out<< std::endl;   

    std::cout<<"内边的比例：" << (1.0* min_innerEdge) / max_innerEdge << std::endl;   

    //打印 ans
    std::cout << "特征数量："<< g.feat_sum.size()<<std::endl;
    std::cout << "Imbalance values: ";

    for (double value : ans) 
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return ans;


}


//检查无误  可以优化为 O(m)
std::tuple<int, std::vector<std::unordered_map<std::string, std::vector<std::pair<int, int>>>>, std::unordered_map<int, std::string>>
calculateGain(const Graph& g, const std::vector<std::vector<int>>& cluster_node_list, std::unordered_map<int, std::vector<int>>& node2neighborsPart, int npart)
{
    // int swap_real_num = 0;
    std::vector<std::unordered_map<std::string, std::vector<std::pair<int, int>>>> dict_list;
    int count_candidate_node = 0;
    std::unordered_map<int, std::string> nodeToKey;

    for (size_t i = 0; i < cluster_node_list.size(); ++i)
    {
        const std::vector<int>& cluster_nodes = cluster_node_list[i];
        std::unordered_map<std::string, std::vector<std::pair<int, int>>> move_dict;
        for (const int& node : cluster_nodes)
        {
            std::vector<int> part_count(npart, 0);
            const std::vector<int>& neighbors = g.edges[node];  // 取出node 的邻居
            for (int neighbor : neighbors)
            {
                int partID = g.node_id_to_partID.at(neighbor);
                ++part_count[partID];
            }
            node2neighborsPart[node] = part_count;

            auto max_iter = std::max_element(part_count.begin(), part_count.end());
            int max_value = *max_iter;
            int max_index = std::distance(part_count.begin(), max_iter);
            int selfPart_value = part_count[g.node_id_to_partID.at(node)];
            int gain = max_value - selfPart_value;
            if (gain > 0) {
                std::string key_str = std::to_string(g.node_id_to_partID.at(node)) + "->" + std::to_string(max_index);
                move_dict[key_str].push_back(std::make_pair(node, gain));
                nodeToKey[node] = key_str;
                ++count_candidate_node;
            }
            else{
                nodeToKey[node] = "";  //gain 小于0 ,node 没有移动方向
            }
        }
        dict_list.push_back(move_dict);
    }

    return std::make_tuple(count_candidate_node, dict_list, nodeToKey);
}



void getIndex(const std::string& key, int& index_1, int& index_2)
{   //std::cout<<" key : "<<key<<std::endl;
    std::stringstream ss(key);
    std::string item;
    std::vector<std::string> parts;
    while (std::getline(ss, item, '>'))
    {
        size_t arrow_pos = item.find('-');
        if (arrow_pos != std::string::npos)
        {
            item.erase(arrow_pos, 1);
        }
        parts.push_back(item);
    }

    // Convert string to int
    index_1 = std::stoi(parts[0]);
    index_2 = std::stoi(parts[1]);
}




//检查无误  省空间版本
std::tuple<int, std::unordered_map<int, std::pair<int,int>>,std::unordered_map<int, int>>
calculateGain(const Graph& g, std::unordered_map<int, std::vector<int>>& node2neighborsPart, const int& npart)
{
    // int swap_real_num = 0;
    // std::vector<std::unordered_map<std::string, std::vector<std::pair<int, int>>>> dict_list;
    int count_candidate_node = 0;
    std::unordered_map<int, std::pair<int,int> > nodeToKey; //只存储gain > 0的点
    std::unordered_map<int, int> nodeToGain; //只存储gain 大于0的点
    for (int i  = 0 ; i< g.n ; i++)
    {

        //const std::vector<int>& cluster_nodes = cluster_node_list[i];
        //std::unordered_map<std::string, std::vector<std::pair<int, int>>> move_dict;
        // for (const int& node : cluster_nodes)
        // {
        int node  = i ;
        std::vector<int> part_count(npart, 0);
        const std::vector<int>& neighbors = g.edges[node];  // 取出node 的邻居
        for (int neighbor : neighbors)
        {
            int partID = g.node_id_to_partID.at(neighbor);
            ++part_count[partID];
        }
        node2neighborsPart[node] = part_count;

        auto max_iter = std::max_element(part_count.begin(), part_count.end());
        int max_value = *max_iter;
        int max_index = std::distance(part_count.begin(), max_iter);
        int selfPart_value = part_count[g.node_id_to_partID.at(node)];
        int gain = max_value - selfPart_value;
        if (gain > 0)
        {
            //std::string key_str = std::to_string(g.node_id_to_partID.at(node)) + "->" + std::to_string(max_index);
            // move_dict[key_str].push_back(std::make_pair(node, gain));
            nodeToKey[node] = std::make_pair(g.node_id_to_partID.at(node), max_index );
            nodeToGain[node] = gain ;
            ++count_candidate_node;
        }
        // else
        // {
        //     nodeToKey[node] = "";  //gain 小于0 ,node 没有移动方向
        // }
        // }
        // dict_list.push_back(move_dict);
        part_count.clear();
        part_count.shrink_to_fit();
    }

    return std::make_tuple(count_candidate_node, nodeToKey, nodeToGain);
}




FiveTuple IfSatisfiedContrainst(const Graph & g , const int& node, const int & index_1, const int & index_2, const std::unordered_map<int, std::vector<int>> & node2neighborsPart, std::vector<int> nodeNumListTmp, std::vector<std::vector<double>> featuresListTmp, std::vector<int> innerEdgeListTmp, std::vector<int> outerEdgeListTmp, double imbalance1, double imbalance2, double imbalance3, double imbalance4, int count_optimal, const std::vector<double>& vector_optimal) //O(1)
{


    bool ans = true;
    int index_1_partNum = node2neighborsPart.at(node)[index_1];
    int index_2_partNum = node2neighborsPart.at(node)[index_2];

    //计算内边
    innerEdgeListTmp[index_1] -= index_1_partNum;
    innerEdgeListTmp[index_2] += index_2_partNum;
    double inerEdgeImbal = static_cast<double>(*std::min_element(innerEdgeListTmp.begin(), innerEdgeListTmp.end())) / *std::max_element(innerEdgeListTmp.begin(), innerEdgeListTmp.end());

    if (inerEdgeImbal < imbalance3)
    {
        //count_ban_innerEdge++;
        ans = false;
    }

    //计算外边line 124line 124
    int index1_out_orig = 0;
    index1_out_orig += (g.edges.at(node).size() - index_1_partNum);

    int index2_out_orig = 0;
    index2_out_orig += index_2_partNum;

    int index1_out_swap = 0;
    index1_out_swap += index_1_partNum;

    int index2_out_swap = 0;
    index2_out_swap += (g.edges.at(node).size() - index_2_partNum);

    outerEdgeListTmp[index_1] += (index1_out_swap - index1_out_orig);
    outerEdgeListTmp[index_2] += (index2_out_swap - index2_out_orig);

    double outerEdgeImbal = static_cast<double>(*std::min_element(outerEdgeListTmp.begin(), outerEdgeListTmp.end())) / *std::max_element(outerEdgeListTmp.begin(), outerEdgeListTmp.end());

    if (outerEdgeImbal < imbalance4)
    {
        //count_ban_outerEdge++;
        ans = false;
    }

    // 计算点数的变化
    nodeNumListTmp[index_1] -= 1;
    nodeNumListTmp[index_2] += 1;
    double nodeNumImbal = 0.0;
    for (int num : nodeNumListTmp) {
        nodeNumImbal = std::max(nodeNumImbal, static_cast<double>(std::abs(num - count_optimal)));
    }
    nodeNumImbal /= count_optimal;

    if (nodeNumImbal > imbalance1) {
        //count_ban_nodeNum++;
        ans = false;
    }

    //计算特征的变化
    std::vector<double> feat = g.nodes.at(node);

    for (size_t i = 0; i < featuresListTmp[index_1].size(); ++i) {
        featuresListTmp[index_1][i] -= feat[i];
    }

    for (size_t i = 0; i < featuresListTmp[index_2].size(); ++i) {
        featuresListTmp[index_2][i] += feat[i];
    }

    std::vector<double> max_abs_diff_divided_2;

    for (size_t i = 0; i < vector_optimal.size(); ++i) {
        max_abs_diff_divided_2.push_back(std::max(std::abs(featuresListTmp[index_1][i] - vector_optimal[i]),
                                                  std::abs(featuresListTmp[index_2][i] - vector_optimal[i])) / vector_optimal[i]);
    }


    double exceed_threshold = *std::max_element(max_abs_diff_divided_2.begin(), max_abs_diff_divided_2.end());

    if (exceed_threshold > imbalance2)
    {
        ans = false;
    }

    int index_1_out_value = index1_out_swap - index1_out_orig;
    int index_2_out_value = index2_out_swap - index2_out_orig;
    int index_1_inner_value = (-index_1_partNum);
    int index_2_inner_value = index_2_partNum;

    return std::make_tuple(ans, index_1_out_value, index_2_out_value, index_1_inner_value, index_2_inner_value);
    
}




void saveResult(std::vector<Cluster>& result, const Graph& g, const int& feacuter_size, const int& k)   //存储result 信息
{

    result.clear();
    for (int i = 0; i < k; ++i) 
    {
        Cluster cluster(0,  std::vector<double>(feacuter_size, 0.0), 0, 0, 0);
        result.push_back(cluster);
    }
    for (const auto& node_data : g.node_id_to_partID)   //这里之前是 g.nodes
    {
        int node = node_data.first;
        // if (g.nodes.find(node) == g.nodes.end()) 
        // {
        //        std::cout << "Key " << node<< " does not exist in the unordered_map g.nodes." << std::endl;
        //        continue;
        // } 


        const std::vector<double>& feat = g.nodes[node];
        int j =  g.node_id_to_partID.at(node);

        int new_count = result[j].count + 1;
        for (size_t i = 0; i < result[j].vector.size(); ++i) {
            result[j].vector[i] = (result[j].vector[i] * result[j].count + feat[i]) / new_count;
        }
        result[j].count = new_count;
        result[j].point[node] =  1 ;
        result[j].nodeSet.insert(node);
        int n_neighbors =  0;
        // if(g.edges.find(node) != g.edges.end())
        // {
        //     n_neighbors = g.edges.at(node).size();
        // }
        // else{
        //     g.edges[node]  =  std::unordered_set<int>();
        // }
        n_neighbors = g.edges[node].size();
        result[j].edgeSum += n_neighbors;
    
    }
    
}



//保存多个移动动作  node--gain --- pair   先小数据集测效果, 后面再优化上大数据集    unordered_map<int, unordered_map< int, pair<int,int>> >
std::tuple<int, std::unordered_map<int, std::pair<int,int>>,std::unordered_map<int, int>>
calculateGainMoreSearchSpace(const Graph& g, std::unordered_map<int, std::vector<int>>& node2neighborsPart, const int& npart)
{
    // int swap_real_num = 0;
    // std::vector<std::unordered_map<std::string, std::vector<std::pair<int, int>>>> dict_list;
    int count_candidate_node = 0;
    std::unordered_map<int, std::pair<int,int> > nodeToKey; //只存储gain > 0的点
    std::unordered_map<int, int> nodeToGain; //只存储gain 大于0的点
    for (int i  =0 ; i< g.n ; i++)
    {

        //const std::vector<int>& cluster_nodes = cluster_node_list[i];
        //std::unordered_map<std::string, std::vector<std::pair<int, int>>> move_dict;
        // for (const int& node : cluster_nodes)
        // {
        int node  =  i ;
        std::vector<int> part_count(npart, 0);
        const std::vector<int>& neighbors = g.edges[node];  // 取出node 的邻居
        for (int neighbor : neighbors)
        {
            int partID = g.node_id_to_partID.at(neighbor);
            ++part_count[partID];
        }
        node2neighborsPart[node] = part_count;

        auto max_iter = std::max_element(part_count.begin(), part_count.end());
        int max_value = *max_iter;
        int max_index = std::distance(part_count.begin(), max_iter);
        int selfPart_value = part_count[g.node_id_to_partID.at(node)];
        int gain = max_value - selfPart_value;
        if (gain > 0)
        {
            //std::string key_str = std::to_string(g.node_id_to_partID.at(node)) + "->" + std::to_string(max_index);
            // move_dict[key_str].push_back(std::make_pair(node, gain));
            nodeToKey[node] = std::make_pair(g.node_id_to_partID.at(node), max_index );
            nodeToGain[node] = gain ;
            ++count_candidate_node;
        }
        // else
        // {
        //     nodeToKey[node] = "";  //gain 小于0 ,node 没有移动方向
        // }
        // }
        // dict_list.push_back(move_dict);
        part_count.clear();
        part_count.shrink_to_fit();
    }

    return std::make_tuple(count_candidate_node, nodeToKey, nodeToGain);
}










void saveResultForIndicator(std::vector<Cluster>& result, Graph& g, const int& feacuter_size, const int& k)   //存储result 信息
{

    result.clear();
    for (int i = 0; i < k; ++i) 
    {
        Cluster cluster(0,  std::vector<double>(feacuter_size, 0.0), 0, 0, 0);
        result.push_back(cluster);
    }
    for (const auto& node_data : g.node_id_to_partID)   //这里之前是 g.nodes
    {
        int node = node_data.first;
        // if (g.nodes.find(node) == g.nodes.end()) 
        // {
        //        std::cout << "Key " << node<< " does not exist in the unordered_map g.nodes." << std::endl;
        //        continue;
        // } 


        const std::vector<double>& feat = g.nodes[node];
        int j =  g.node_id_to_partID.at(node);

        int new_count = result[j].count + 1;
        for (size_t i = 0; i < result[j].vector.size(); ++i) {
            result[j].vector[i] = (result[j].vector[i] * result[j].count + feat[i]) / new_count;
        }
        result[j].count = new_count;
        result[j].point[node] =  1 ;
        result[j].nodeSet.insert(node);
        int n_neighbors =  0;
        // if(g.edges[node] != g.edges.end())
        // {
        //     n_neighbors = g.edges.at(node).size();
        // }
        // else{
        //     g.edges[node]  =  std::unordered_set<int>();
        // }
        
        n_neighbors = g.edges[node].size();
        result[j].edgeSum += n_neighbors;
    
    
    
    }

    
    
}

