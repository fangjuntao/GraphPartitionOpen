//BIHM(节点按照degree 大小从大到小进行排序)
#include "graph.h"
#include "cluster.h"
#include "evaluation.h"
#include <set>
#include<unordered_set>
#include<unordered_map>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib> // for  srand
#include "omp.h"
#include <limits>
#include <iostream>
#include <fstream>
#include "GainList.h"
#include "kmeans/kmeans.hpp"
using namespace std;
//vector 遍历速度远大于 set 
//没必要存string 存pair 也可以
//feature 不要一开始就存avrage 会有精度损失
//constraints node set 实现真实的O(1)修改和删除，插入:  通过链表去维护一个队列，存好所在的位置即可
//可以优化的地方标记, 可以优化


using FourTuple = std::tuple<int, int, int, int>;
using FiveTuple = std::tuple<bool, int, int, int, int>;


// Function to normalize the features of nodes
vector<vector<double>> normalize_features(const std::vector<std::vector<double>>& nodes_orign) 
{

    vector<vector<double>> nodes  = nodes_orign ;
    int nobs = nodes.size();  // 数据点数量
    int ndim = nodes[0].size();  // 数据维度

    // Find min and max for each dimension
    std::vector<double> min_vals(ndim, std::numeric_limits<double>::max());
    std::vector<double> max_vals(ndim, std::numeric_limits<double>::lowest());

    for (const auto& node : nodes) 
    {
        for (int j = 0; j < ndim; ++j) 
        {
            min_vals[j] = std::min(min_vals[j], node[j]);
            max_vals[j] = std::max(max_vals[j], node[j]);
        }
    }

    // Normalize the nodes
    for (auto& node : nodes) {
        for (int j = 0; j < ndim; ++j) 
        {
            if (max_vals[j] != min_vals[j]) 
            {
                node[j] = (node[j] - min_vals[j]) / (max_vals[j] - min_vals[j]);
            } else 
            {
                node[j] = 0.0; // if all values are the same, set to 0
            }
        }
    }



    return nodes ;

}





vector<int>  kmeansGraph(const Graph& g , int k)
{


/*     // 假设我们有以下数据矩阵（示例数据）
    std::vector<std::vector<double>> nodes = {
        {1.0, 2.0},
        {1.5, 1.8},
        {5.0, 8.0},
        {8.0, 8.0},
        {1.0, 0.6},
        {9.0, 11.0},
        {8.0, 2.0},
        {10.0, 2.0},
        {9.0, 3.0}
    };
 */
    // 归一化节点特征
    
    vector<vector<double>> nodes  =  normalize_features(g.nodes);

    int nobs = nodes.size();  // 数据点数量
    int ndim = nodes[0].size();  // 数据维度
    int ncenters = k;  // 簇数量

    // 将 nodes 转换为一维向量，以适应 SimpleMatrix 的要求
    std::vector<double> matrix;
    for (const auto& node : nodes) 
    {
        matrix.insert(matrix.end(), node.begin(), node.end());
    }

    // 将矩阵封装在 SimpleMatrix 中
    kmeans::SimpleMatrix<double, int, int> kmat(ndim, nobs, matrix.data());

    // 计算 k-means 聚类
    auto res = kmeans::compute(
        kmat,
        kmeans::InitializeKmeanspp<kmeans::SimpleMatrix<double, int, int>, int, double>(),  // 使用 kmeans++ 初始化
        kmeans::RefineLloyd<kmeans::SimpleMatrix<double, int, int>, int, double>(),  // 使用 Lloyd 算法精炼
        ncenters
    );

    vector<int> ans  = res.clusters; //聚类结果


    
/*     // 输出聚类结果
    std::cout << "Cluster centers:\n";
    for (int i = 0; i < ncenters; ++i) 
    {
        for (int j = 0; j < ndim; ++j) 
        {
            std::cout << res.centers[j + i * ndim] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nCluster assignments:\n";
    for (int i = 0; i < nobs; ++i) 
    {
        std::cout << res.clusters[i] << " ";
    }
    std::cout << "\n";

    // 计算 WCSS
    std::vector<double> wcss(ncenters);
    kmeans::compute_wcss
    (
        kmat, 
        ncenters, 
        res.centers.data(), 
        res.clusters.data(), 
        wcss.data()
    );

    std::cout << "\nWCSS:\n";
    for (double val : wcss) 
    {
        std::cout << val << " ";
    }
    std::cout << "\n"; */


    return ans;

}


void BalancedInitialization(Graph &g, const vector<int> & kmeansResult, int  k , int p ) //均匀初始化
{
    vector<int> clusterCount(k,0);
    unordered_map<int,bool> vertexState;  //true 代表 assigned      flase 代表的是unassigned
    for(int i = 0 ; i < g.n ; i++)
    {
        clusterCount[kmeansResult[i]]++ ;
        vertexState[i]   =  false;
    }

    vector<int>partEdge(p, 0);
    vector<int>partCut(p, 0);
    vector<vector<int>> partVertex;
    vector<int> clusterCountTmp(k,0);
    for(int i = 0 ;  i < p ; i++)
    {
        partVertex.push_back(clusterCountTmp);
    }
    int totalEdge =  0 ;
    int totalCut  =  0 ;

    // 创建一个包含节点索引及其邻居数量的 vector
    std::vector<std::pair<int, int>> nodeNeighborCount;
    for (int i = 0; i < g.edges.size(); ++i) {
        nodeNeighborCount.push_back({i, g.edges[i].size()});
    }

    // 按邻居数量从大到小排序
    std::sort(nodeNeighborCount.begin(), nodeNeighborCount.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) 
            {
                return a.second > b.second;
            });

    // nodeNeighborCount 现在按邻居数量从大到小排序
    // 可以打印或使用排序后的节点顺序
    for (const auto& node_pair : nodeNeighborCount)
    {
        int node =  node_pair.first;
        int visitNeighbor = 0 ;
        vector<int> neiPart(p, 0);
        vector<int> score(p, 0);
        for(const auto& u : g.edges[node])
        {
            if( vertexState[u])
            {
                neiPart[g.node_id_to_partID[u]]++;
                visitNeighbor++;
            }
        }
        int degreeOut, edge,cut, partSize, totalNode;


        // 计算score
        for (int i = 0; i < p; ++i)
        {
            degreeOut =  visitNeighbor - neiPart[i];
            edge  =  partEdge[i] + neiPart[i];
            cut  = partCut[i]+ degreeOut ;
            partSize  = partVertex[i][kmeansResult[node]]+1;
            totalNode  =  clusterCount[kmeansResult[node]];
            score[i]  =  partSize;

            if (visitNeighbor != 0) 
            {
                score[i] += degreeOut / visitNeighbor;
            }

            if (totalEdge != 0 && partEdge[i] > (totalEdge / p))
            {
                if (neiPart[i] > 0) 
                {
                    score[i] += 1;
                }
            }

            if ( totalCut != 0 && partCut[i] > (totalCut / p)) {
                if (degreeOut > 0) 
                {
                    score[i] += 1;
                }
            }
        }

    
        // 使用 std::min_element 找到最小值的迭代器
        auto minElementIt  = std::min_element(score.begin(), score.end());

        // 计算下标
        int selectedPartID = std::distance(score.begin(), minElementIt);

        g.node_id_to_partID[node]   =  selectedPartID ;


        /* Update 相关变量*/

        vertexState[node] =  true ;
        partEdge[selectedPartID ]   =  partEdge[selectedPartID]+neiPart[selectedPartID];
        partCut[selectedPartID]     = partCut[selectedPartID]+visitNeighbor-neiPart[selectedPartID];
        partVertex[selectedPartID][kmeansResult[node]]++;
        totalEdge += neiPart[selectedPartID];
        totalCut += (visitNeighbor-neiPart[selectedPartID]) ;

        for(int i =  0 ; i< neiPart.size() ; i++)
        {

            if(i != selectedPartID)
            {
                partCut[i] += neiPart[i];
            }
        }

    }





}


//最好也再检查一下
void moveNodesBin(Graph &g,  unordered_map<int,   pair<int,int> >& nodeToKey , unordered_map<int,  int>& nodeToGain, vector<int>& nodeNumList, vector<vector<double>>& featuresList, vector<int>& innerEdgeList, vector<int>& outerEdgeList, unordered_map<int,  vector<int>> &node2neighborsPart, int k,  int args_swap_sum, double imbalance1, double imbalance2, double imbalance3, double imbalance4, int count_optimal, const vector<double>& vector_optimal ) //省空间版本，加速
{
  
    unordered_map<int, FourTuple> node_FourTuple ;  // （node, value1,value2,value3,value4)
    //注意以下二者的nodeID gain的对应关系需要维护好

   // unordered_map<int,  int> nodeToGain ;   //nodeID -gain
    unordered_map<int, bool> nodeToConstraints ;  // nodeID--bool  true 代表满足constrainsts  false 代表不满足contrainsts，只存 true 节点

    vector<pair<int,int>> insertVector;  //index, value
    vector<pair<int,int>> deleteVector ;  //index,value


    //计算graph的最大degree

    GainList gainList(10000);  //初始化大小
    
    int swap_real_num = 0; //记录真实移动的次数

    // vector<int> max_gain_index ;
    int max_gain_sum = 0 ; //记录移动的minimize total cut 的数量

    int max_gain               = 0 ;
    int count_ban_innerEdge    = 0 ;
    int count_ban_outerEdge    = 0 ;
    int count_ban_nodeNum      = 0 ;
    int count_ban_features     = 0 ;
    int count_candidate_node   = 0 ;

    int count_modifiNodes =  0;
    //初始化变量
    int index_1_out_value      = -1 ; 
    int index_2_out_value      = -1 ;
    int index_1_inner_value    = -1 ;
    int index_2_inner_value    = -1 ; 

    for(int i = 0 ; i< g.n ; i++)   //O(n)    unodered_map的遍历顺序无法保证
    {
        int  node        = i;
        int  gain        =  0 ; 
        if( nodeToGain.find(node) != nodeToGain.end() )
        {
            gain = nodeToGain.at(node);  //gain 必不为0
        }
        else
        {
            continue;
        }


        // auto it = cluster_dict.begin();
        // std::advance(it, i);
        // auto kv =  it ;
        // const  string &key   = kv->first;  //string 移动的方向
        // const  auto   &value = kv->second;  // vector 里面存的是（node, gain)


        int index_1 = nodeToKey.at(node).first ;
        
        int index_2 = nodeToKey.at(node).second;
       // getIndex(key,index_1,index_2) ;   // Convert string to int  pair<int, int>
        
        
        
        // auto node_gain   = value[i];

        
        bool  constraints ; 
        tie(constraints, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
            IfSatisfiedContrainst(g, node, index_1,  index_2, node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

        node_FourTuple[node]  =  make_tuple(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value);  
        //bool constraints  =   IfSatisfiedContrainst(node, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;
        

       // nodeToGain[node]   =   gain  ;  //默认为0，表示不在set 里面

        if( !constraints )  //不满足contraints
        {
            // nodeToGain[node]  =  -1  ;  //不满足constraints的gain 初始化为-1 ; 
            continue;
        }
        else //满足contraints
        {
            nodeToConstraints[node]  =  constraints ;   //只存满足contrainst的节点
            // nodeToGain[node]  =  gain ; //记录节点的gain
            if (gain > 0 )  //gain 一定大于0
            {
                // auto it = constraintsNodeSet[gain].find(node);
                // if (it != constraintsNodeSet[gain].end()) 
                // {
                //     std::cout << "元素 " << node << " 存在于vector中。" << std::endl;
                //     cout<<" 寄了寄了 line 136 "<<endl;
                // }
                // if(node ==458752 )
                // {
                //     cout<<"line 130 : ";
                //     cout<<"nodeToGain[458752]: "<<nodeToGain[node]<<endl;
                //     cout<<"nodeToKey[458752]:"<<nodeToKey[458752] <<endl;
                // }
                //constraintsNodeSet[gain].insert(node);

                gainList.insert(Node(node,gain));
                max_gain   = max(max_gain, gain);
                ++count_candidate_node ;
            }
        
        }


    
    
    
    }

    // cout<< "max_gain: "<< max_gain<<endl;
    
    

    int findMaxNum  =  0 ;
    while(1)         //开始进行move
    {
        int swap_real_numCopy  =  swap_real_num;

        for(int i = max_gain; i > 0; --i )    //存在gain > 0 且满足constraints的点
        {
            deleteVector.clear();
            insertVector.clear();
           // for(const int& node_1: constraintsNodeSet[i]  )  //当您在遍历 set 时尝试修改它（如删除元素）时，可能会导致未定义的行为。在这种情况下，您应该使用另一种方法来删除元素，例如使用迭代器。
            for(const Node& node_1 : gainList[i])
            {            //auto it = constraintsNodeSet[i].begin(); it != constraintsNodeSet[i].end();     
                // Node max_gain_node  = node_heap.elements[0];

                int  gain    =   i ; 
                
                bool  constraintsThis;


                int index_1 =  nodeToKey[node_1.id].first;
                int index_2 =  nodeToKey[node_1.id].second;
                // getIndex(key, index_1, index_2);
                tie(constraintsThis, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
                    IfSatisfiedContrainst(g, node_1.id, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

                if( constraintsThis == false)
                {
                    //删除元素

                    --count_candidate_node  ;
                    ++findMaxNum ;
                    //constraintsNodeSet[i].erase(node_1); 
                    nodeToConstraints.erase(node_1.id);
                    // nodeToGain.erase(node_1.id);
                    // nodeToKey.erase(node_1.id);
                    
                    deleteVector.push_back(make_pair(i,node_1.id)) ; 
                    continue;  

                }
                else
                {

                    max_gain  = gain ;  


                    ++swap_real_num ;
                    // 访问一个映射
                    int index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value;
                    tie(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = node_FourTuple.at(node_1.id);


                    //更新相关信息
                    g.node_id_to_partID[node_1.id] = index_2;

                    innerEdgeList[index_1] += index_1_inner_value;
                    innerEdgeList[index_2] += index_2_inner_value;
                    outerEdgeList[index_1] += index_1_out_value;
                    outerEdgeList[index_2] += index_2_out_value;

                    nodeNumList[index_1] -= 1;
                    nodeNumList[index_2] += 1;


                    const vector<double>& feat = g.nodes[node_1.id];

                    for (size_t i = 0; i < featuresList[index_1].size(); ++i) {
                        featuresList[index_1][i] -= feat[i];
                    }

                    for (size_t i = 0; i < featuresList[index_2].size(); ++i) {
                        featuresList[index_2][i] += feat[i];
                    }

                    // //更新node_1 本身的信息
                    // vector< pair<int, int>> new_vector;

                    // for (const auto &pair : dict_list[g.node_id_to_clusterID.at(node_1.id) ].at(nodeToKey[node_1.id])) 
                    // {
                    //     if (pair.first != node_1.id) 
                    //     {
                    //         new_vector.push_back(pair);
                    //     }
                    // }
                    // dict_list[g.node_id_to_clusterID.at(node_1.id) ].at(nodeToKey[node_1.id]) = new_vector;
                    // if (new_vector.size() == 0)
                    // {
                    //     dict_list[g.node_id_to_clusterID.at(node_1.id)].erase(nodeToKey[node_1.id]);
                    // }
                    int returnValue1  = nodeToKey.erase(node_1.id) ;
                    int returnValue2 =  nodeToGain.erase(node_1.id);
                    if(returnValue1 == 0  || returnValue2 == 0 )
                    {
                        cout<<" wrong in main line 245***************"<<endl;
                    }

                    //在set 里面删除node_1
                    //constraintsNodeSet[gain].erase(node_1);
                   // constraintsNodeSet[0].insert(node_1);  意义不大，还增加负担

                    deleteVector.push_back(make_pair(gain, node_1.id));
                    nodeToGain.erase(node_1.id);
                    nodeToKey.erase(node_1.id);
                    nodeToConstraints.erase(node_1.id);

                    // //更新邻居的信息
                    for(const auto & neighbor: g.edges.at(node_1.id))
                    {
                        // if(neighbor ==458752 )
                        // {
                        //     cout<<"line 370 : ";
                        //     cout<<"nodeToGain[458752]: "<<nodeToGain[neighbor]<<endl;
                        //     cout<<"nodeToKey[458752]:"<<nodeToKey[458752] <<endl;
                        // }
                        int node1_original_partID  = index_1 ;
                        int node1_swap_to_partID   =  index_2 ; 
                        // int index = g.node_id_to_clusterID.at(neighbor);
                        node2neighborsPart[neighbor][node1_original_partID] = max(node2neighborsPart.at(neighbor)[node1_original_partID]-1,0);
                        node2neighborsPart[neighbor][node1_swap_to_partID] += 1;

                        
                        auto max_iter = max_element(node2neighborsPart.at(neighbor).begin(), node2neighborsPart.at(neighbor).end());
                        int max_value = *max_iter;
                        int max_index = distance(node2neighborsPart.at(neighbor).begin(), max_iter);
                        

                        int  selfPart_value = node2neighborsPart.at(neighbor)[g.node_id_to_partID.at(neighbor)];
                        int  gain =  max_value - selfPart_value ;
                        int index_1_tmp  =  g.node_id_to_partID.at(neighbor) ;



                        int index_2_tmp  = max_index ;
                        bool  constraints ; 
                        tie(constraints, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
                        IfSatisfiedContrainst(g, neighbor, index_1_tmp,  index_2_tmp,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

                        node_FourTuple[neighbor]     =  make_tuple(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value);  
                    
                        //string  move_str =    nodeToKey[neighbor]  ;     //move_str 如果为空说明 之前这个点不在move_dict 里面
                       // string  key_str  =    to_string(g.node_id_to_partID.at(neighbor)) + "->" +    to_string(max_index);
                        // if( gain > 0) 
                        // {
                        //     // getIndex(key_str, index_1_tmp, index_2_tmp);
                        //         // if( index_1_tmp == index_2_tmp)  //二者可能相等  gain==0
                        //         // {
                        //         // }
                        // }


                        if ( nodeToGain.find(neighbor) != nodeToGain.end() && nodeToConstraints.find(neighbor) !=nodeToConstraints.end() )  //之前在vector里面 
                        {
                            

                            // auto it  = constraintsNodeSet[nodeToGain[neighbor]].find(neighbor);
          
                            // if (it == constraintsNodeSet[nodeToGain[neighbor]].end()) 
                            // {

                            //         cout<<" cuolecuole 338 "<<endl;

                            // }
                            //判断现在是否满足从constrains

                            // 1. gain发生改变
                            if(nodeToGain[neighbor] != gain)
                            {


                                //现在的gain 为 0  或者现在不满足contraints
                                if(gain == 0 ||constraints== false ) //在原来的vector里面删除它
                                {
                                // auto it_map = node_set.find(neighbor);  //注意需要是nodeID 唯一才可以，就是里面的nodeID只能在set出现一次
                                    //node_set.erase({neighbor,nodeToGain[neighbor]});  
 
                                    //node_heap.Delete(neighbor); //删除元素
                                    // if(neighbor ==458752 )
                                    // {
                                    //     cout<<"line 370 : ";
                                    //     cout<<"nodeToGain[458752]: "<<nodeToGain[neighbor]<<endl;
                                    //     cout<<"nodeToKey[458752]:"<<nodeToKey[458752] <<endl;
                                    // }
                                    //constraintsNodeSet[nodeToGain[neighbor]].erase(neighbor);

                                    deleteVector.push_back(make_pair(nodeToGain[neighbor],neighbor));

                                    
                                    --count_candidate_node ;

                                
                                }
                                else  //gain 大于 0 且满足constraints
                                {
                                    //1. 将原来删除 
                                    // auto it_map = node_set.find(neighbor);

                                   // constraintsNodeSet[nodeToGain[neighbor]].erase(neighbor);

                                    deleteVector.push_back(make_pair(nodeToGain[neighbor],neighbor));

                                    //2. 将现在插入

                                // insert_node(neighbor, gain );  //插入节点,包括同步存储 node - gain

                                 //   node_heap.Push(Node(neighbor, gain));
                                 //   constraintsNodeSet[gain].insert(neighbor);
                                    insertVector.push_back(make_pair(gain, neighbor));
                                    max_gain  =  max(gain, max_gain);
                            
                                
                                }

                            }
                            else  //gain 不变
                            {
            
                                //contrainst 发生变化
                                if( constraints != nodeToConstraints[neighbor])
                                {
                                    // if( !nodeToConstraints[neighbor])
                                    // {
                                    //     cout<<"wrong wrong line 620 "<<endl;
                                    // }
                                    // auto it_map = node_set.find(neighbor);
                                    //node_set.erase({neighbor,nodeToGain[neighbor]});

                                    //node_heap.Delete(neighbor);
                                   // constraintsNodeSet[nodeToGain[neighbor]].erase(neighbor);

                                    deleteVector.push_back(make_pair(nodeToGain[neighbor], neighbor) );
                                    --count_candidate_node ;

                                }


                            }

                        }
                        else  //不在heap 里面
                        {

                            //现在需要进去set
                            if( gain > 0 && constraints)
                            {


                                
                                //constraintsNodeSet[gain].insert(neighbor);
                                insertVector.push_back(make_pair(gain,neighbor ));
                                max_gain  =  max(gain, max_gain);
                                count_candidate_node ++;
                            }

                            //现在依旧不需要进入set

                        

                        }

                    
                        if(gain > 0 ) 
                        {
                            nodeToGain[neighbor]  =  gain ;
                            nodeToKey[neighbor]   =   make_pair(g.node_id_to_partID.at(neighbor),max_index) ;
                        }
                        else{
                            
                            nodeToGain.erase(neighbor);
                            nodeToKey.erase(neighbor);

                        }
                        if(constraints == true)
                        {
                            nodeToConstraints[neighbor]  = constraints ;
                        }
                        else
                        {
                            nodeToConstraints.erase(neighbor); 
                        }
                    
                    }

                    //更新contraints 可能发生变化的节点    X -> index_1  and  index_1 -> X  and X-> index_2 index_2 -> X

                    max_gain_sum += max_gain;


                    if (swap_real_num % 100 == 0)
                    {
                        
                        cout<<"swap_real_num: "<<swap_real_num<<endl;
                        
                        // 使用    chrono::system_clock 获取当前时间
                        auto current_time =    chrono::system_clock::now();
                        
                        // 将当前时间转换为    time_t 类型

                        time_t time_t_current_time =    chrono::system_clock::to_time_t(current_time);

                        // 使用    ctime 将    time_t 类型的当前时间转换为字符串
                        string time_str =    ctime(&time_t_current_time);
                        
                        // 打印当前时间字符串

                        cout <<"max_gain:"<<max_gain<<endl;
                        //并行的话会不准
                        cout<<" max_gain_sum :"<<max_gain_sum<<endl;
                        cout<<"count_ban_innerEdge:"<<count_ban_innerEdge<<endl;
                        cout<<"count_ban_outerEdge:"<<count_ban_outerEdge<<endl;
                        cout<<"count_ban_feauture:"<< count_ban_features<<endl;
                        cout<<"count_ban_nodeNum: "<<count_ban_nodeNum<<endl;
                        cout<<"count_candidate_node: "<<count_candidate_node<<endl;
                        cout<<" count_modifiNodes: "<<count_modifiNodes<<endl;
                        cout<<" 100 轮总共查找最大值的次数："<<findMaxNum<<endl;
                        cout << "Current time: " << time_str << endl;
                        findMaxNum  = 0 ;
                        count_modifiNodes = 0 ;


                    }


                    i =  0 ;
                }

                // node_heap.Pop();   //删除第一个点

                break;
            }

            count_modifiNodes += deleteVector.size();
            for(auto &pair: deleteVector)
            {

                int index =  pair.first;
                int node  =  pair.second;
                // if (constraintsNodeSet[index].find(node) != constraintsNodeSet[index].end()) 
                // {
                    // 删除元素
                    //constraintsNodeSet[index].erase(node);
                gainList.remove(node);
               // } 
                // else 
                // {
                //      cout << "Value: " << node << " does not exist in the set.\n";
                // }

            }
            count_modifiNodes += insertVector.size();
            for(auto& pair:insertVector)
            {

                int index =  pair.first;
                int node  =  pair.second;

                //constraintsNodeSet[index].insert(node);
                gainList.insert(Node(node,index));
                   
                //}

            }
        
        }

        if(swap_real_numCopy == swap_real_num)
        {
            cout<<" 本轮没有move "<<endl;
            // constraintsNodeSet.clear();
            // constraintsNodeSet.resize(g.n);
            gainList.clear();
            gainList.resize(g.n);
            max_gain =  0;
            count_candidate_node  = 0 ;
            // dict_list.clear();
            //O(n) 进行更新  gainList
                    
            for(int i  =0 ; i< g.n ; i++)   //O(n)    
            {
                int  node        =  i ;
                int  gain        =  0 ; 
                if( nodeToGain.find(node) != nodeToGain.end() )
                {
                    gain = nodeToGain.at(node);  //gain 必不为0
                }
                else
                {
                    continue;
                }

                // auto it = cluster_dict.begin();
                // std::advance(it, i);
                // auto kv =  it ;
                // const  string &key   = kv->first;  //string 移动的方向
                // const  auto   &value = kv->second;  // vector 里面存的是（node, gain)


                int index_1 = nodeToKey.at(node).first ;
                
                int index_2 = nodeToKey.at(node).second;
                // getIndex(key,index_1,index_2) ;   // Convert string to int  pair<int, int>
                
                
                
                // auto node_gain   = value[i];

                
                bool  constraints ; 
                tie(constraints, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
                    IfSatisfiedContrainst(g, node, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

                node_FourTuple[node]  =  make_tuple(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value);  
                //bool constraints  =   IfSatisfiedContrainst(node, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;
                

                // nodeToGain[node]   =   gain  ;  //默认为0，表示不在set 里面

                if( !constraints )  //不满足contraints
                {
                    // nodeToGain[node]  =  -1  ;  //不满足constraints的gain 初始化为-1 ; 
                    continue;
                }
                else //满足contraints
                {
                    nodeToConstraints[node]  =  constraints ;   //只存满足contrainst的节点
                    // nodeToGain[node]  =  gain ; //记录节点的gain
                    if (gain > 0 )  //gain 一定大于0
                    {
                        // auto it = constraintsNodeSet[gain].find(node);
                        // if (it != constraintsNodeSet[gain].end()) 
                        // {
                        //     std::cout << "元素 " << node << " 存在于vector中。" << std::endl;
                        //     cout<<" 寄了寄了 line 136 "<<endl;
                        // }
                        // if(node ==458752 )
                        // {
                        //     cout<<"line 130 : ";
                        //     cout<<"nodeToGain[458752]: "<<nodeToGain[node]<<endl;
                        //     cout<<"nodeToKey[458752]:"<<nodeToKey[458752] <<endl;
                        // }
                        //constraintsNodeSet[gain].insert(node);

                        gainList.insert(Node(node,gain));
                        max_gain   = max(max_gain, gain);
                        ++count_candidate_node ;
                    }
                
                }


            
            
            
            }

            // cout<< "max_gain: "<< max_gain<<endl;
            
            if(max_gain == 0)
            {

                cout<<"无法继续move"<<endl;
                break;
            }

        }


    }
        cout<<" max_gain_sum: "<<max_gain_sum<<endl;

        cout<<"swap_real_num: "<<swap_real_num<<endl;  //实际上move 次数

}







void write_map_to_file(const Graph &g, const std::string& filename) {
    std::ofstream outfile(filename, std::ios::trunc);  // Open the file in truncate mode to clear it before writing

    if (outfile.is_open()) {
        for (const auto& pair : g.node_id_to_partID) {
            int nodeId  = pair.first;
            outfile << g.graph_id_to_orig_id.at(nodeId)  << "," << pair.second << "\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open file";
    }
}


/*


int main() 
{


   string data_input = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Pokec social network/data1,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Pokec social network/data2,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Pokec social network/data3";

    //参数
    //string data_input  = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data1,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data2,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data3";
    //string data_input = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data1,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data2,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data3";
   //string data_input =  "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data1,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data2,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data3";
    // unsigned int seed  =  50;

    int  k  =  10 ;  
    int    feacuter_size = 2  ;   //注意不同的数据集需要手动设置
    double imbalance1    = 0.10  ; //node
    double imbalance2    = 0.10  ; //feature
    double imbalance3    = 0.60; //innerEdge
    double imbalance4    = 0.60;  //outerEdge
    unsigned int seed    =  50; 
    const int    args_move_sum  = 200000;
    const int    args_swap_sum  = 200000;
    // stringstream(argv[2]) >> imbalance1;
    // stringstream(argv[3]) >> imbalance2;
    // stringstream(argv[4]) >> imbalance3;
    // stringstream(argv[5]) >> imbalance4;
    // stringstream(argv[6]) >> seed;
    // stringstream(argv[7]) >> feacuter_size;
    // stringstream(argv[8]) >> k;
    double nodeImbalnce       =  numeric_limits<double>::infinity();  //无穷大
    double featureImbalance   =  numeric_limits<double>::infinity();   //无穷大
    double innerEdgeImbance   =  0.0  ; 
    double outerEdgeImbalance =  0.0  ;
    cout<<"seed: "<<seed<<endl;
    cout<<"data_input: "<<data_input<<endl;
    cout<<"imbalance1:"<<imbalance1<<endl;
    cout<<"imbalance2:"<<imbalance2<<endl;
    cout<<"imbalance3:"<<imbalance3<<endl;
    cout<<"imbalance4:"<<imbalance4<<endl;
    srand(seed);

    // string data_input  =  "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data3/part-000.csv";
    // string data_input  = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data3/part-000.csv";
    
    //string data_input =  "/data/home/juntaofang/jxqy/data1/part-00000.csv,/data/home/juntaofang/jxqy/data2/part-00000.csv,/data/home/juntaofang/jxqy/data3/part-00000.csv"; 
    Graph g = read_graph(data_input);  //可以优化

    // Do something with the graph
    cout<<" n: "<<g.n<<endl;
    cout<<" m: "<<g.m<<endl; //只有适合于文件把边存了的情况

    double count_optimal  = g.n/k ;  //最佳点数
    vector<double> vector_optimal(g.feat_avg.size());
    for (size_t i = 0; i < g.feat_avg.size(); ++i) 
    {
        vector_optimal[i] = g.feat_avg[i] * g.n/k;
    }

    
    vector<Cluster> result;  //save result  可以优化---> 空间

    vector<int> nodeNumList;
    vector<vector<double>> featuresList;
    vector<long long > innerEdgeList ;
    vector<long long> outerEdgeList ; 

    double time_start = (double) clock();

    double time_random_start = (double) clock();
    //随机分配，直到达到要求
    while ( nodeImbalnce > imbalance1 || featureImbalance > imbalance2|| innerEdgeImbance < imbalance3 || outerEdgeImbalance < imbalance4)
    {

        result.clear();
        result.resize(k, Cluster(0,  vector<double>(feacuter_size, 0.0), 0, 0, 0));

        // Iterate over nodes in the graph
        for (const auto& node_data : g.nodes) 
        {
            int node = node_data.first;

            const  vector<double>& feat = node_data.second;

            int part_num = rand() % k; // Randomly assign partition
            int j = part_num;
            g.node_id_to_partID[node] = part_num ; // Store partID in the node data

            int new_count = result[j].count + 1;
            for (size_t i = 0; i < result[j].vector.size(); ++i) {
                result[j].vector[i] = (result[j].vector[i] * result[j].count + feat[i]) / new_count;  //存总数
            }
            result[j].count = new_count;
            result[j].point[node]  =  1 ;  //不为0 代表partition包含这个节点 shanchu
            result[j].nodeSet.insert(node);  
            int n_neighbors = g.edges[node].size();
            result[j].edgeSum += n_neighbors;

        }


        vector<double> ans =   computeObject(result,nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  k );

        nodeImbalnce  =  ans[0];
        double feutureThis = 0.0 ;
        for(int i = 0 ;  i< g.feat_avg.size(); ++i) 
        {
            if(ans[1+i]> feutureThis)  feutureThis =  ans[1+i] ; //存多维特征里面最不均匀的那个特征的imbalance
        }
        featureImbalance    =  feutureThis ; //存多维特征里面最不均匀的那个特征的imbalance
        innerEdgeImbance    =  ans[1+g.feat_avg.size()];
        outerEdgeImbalance  =  ans[2+g.feat_avg.size()]; 

    }


    double time_random_end = (double) clock();
    //k-means先跳过


    //vector<vector<int>> cluster_node_list(1); // 这里由于没有k-means， 所以所有点都是处于同一层   可以优化---> 空间

    size_t i = 0;
    // for (const auto& node_pair : g.nodes) 
    // {
    //     int node = node_pair.first;

    //     cluster_node_list[0].push_back(node);

    //     ++i;
    // }

    //1. 计算gain
    int count_candidate_node;
  //  vector< unordered_map<  string,   vector< pair<int, int>>>> dict_list;
    unordered_map<int,   pair<int,int>> nodeToKey;  //只存gain > 0的点
    unordered_map<int,   int> nodeToGain;   //只存gain > 0的点
    unordered_map<int,  vector<int>> node2neighborsPart;

    double time_gain_start = (double) clock();

    tie(count_candidate_node,  nodeToKey, nodeToGain) = calculateGain( g, node2neighborsPart , k);  //可以优化---> 时间 优化为O(m)
    double time_gain_end = (double) clock();

    cout<<"1125 args_swap_sum: "<<args_swap_sum<<endl;

    //进入move 阶段

    double time_move_start = (double) clock();
    cout<<" 1435 "<<endl ;
    moveNodesBin(g, nodeToKey, nodeToGain,nodeNumList,featuresList,innerEdgeList,outerEdgeList , node2neighborsPart, k,  args_move_sum, imbalance1, imbalance2, imbalance3,  imbalance4, count_optimal, vector_optimal );
   
   
    saveResult(result ,  g, feacuter_size,k ) ;

    vector<double> ans_move =   computeObject(result, nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  k );


    double time_move_end = (double) clock();
    
   
   
    // 进入swap阶段
    // double time_swap_start = (double) clock();
    // swap_pairs(g , nodeToKey ,cluster_node_list,    dict_list,  
    // node2neighborsPart,  nodeNumList, 
    // featuresList, innerEdgeList,  outerEdgeList, 
    // k,  args_swap_sum, imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  );
 
    
    // saveResult(result ,  g, feacuter_size,k ) ;

    // vector<double> ans_swap =   computeObject(result, nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  k );

    // double time_swap_end = (double) clock();

    
    //把结果写入文件系统

    string filename = "output.txt";

    write_map_to_file(g, filename);


    double time_end = (double) clock();
    cout << "Algorithm Time: " << (time_end- time_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Random Ininital Time: " << (time_random_end- time_random_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Compute Gain  Once Time: " << (time_gain_end- time_gain_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Move Time: " << (time_move_end - time_move_start)/CLOCKS_PER_SEC  << "s" << endl;
    //cout << "Swap Time: " << (time_swap_end - time_swap_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "over!!"<<endl;
    return 0;

}



 */




int main(int argc, char* argv[]) 
{

    if(argc < 11) 
    {
        cout << "请提供数据输入路径、imbalance1, imbalance2, imbalance3, imbalance4, seed, feacuter_size, k, p, output 作为命令行参数." << endl;
        return 1;
    }
    string data_input  = argv[1];
    string data_output = argv[10];
    
    //参数
    //string data_input  = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/graphPartition/facebook/data3/part-000.csv";
    //string data_input = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data3/part-000.csv";
    //string data_input =  "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data3/part-000.csv";
    // unsigned int seed  =  50;

    int  k, p ;  
    int    feacuter_size  ;   //注意不同的数据集需要手动设置
    double imbalance1  ; //node
    double imbalance2  ; //feature
    double imbalance3  ; //innerEdge
    double imbalance4  ;  //outerEdge
    unsigned int seed  ; 
    const int    args_move_sum  = 200000;
    const int    args_swap_sum  = 200000;
    stringstream(argv[2]) >> imbalance1;
    stringstream(argv[3]) >> imbalance2;
    stringstream(argv[4]) >> imbalance3;
    stringstream(argv[5]) >> imbalance4;
    stringstream(argv[6]) >> seed;
    stringstream(argv[7]) >> feacuter_size;
    stringstream(argv[8]) >> k;
    stringstream(argv[9]) >> p;
    double nodeImbalnce       =  numeric_limits<double>::infinity();  //无穷大
    double featureImbalance   =  numeric_limits<double>::infinity();   //无穷大
    double innerEdgeImbance   =  0.0  ; 
    double outerEdgeImbalance =  0.0  ;
    cout<<"seed: "<<seed<<endl;
    cout<<"data_input: "<<data_input<<endl;
    cout<<"imbalance1:"<<imbalance1<<endl;
    cout<<"imbalance2:"<<imbalance2<<endl;
    cout<<"imbalance3:"<<imbalance3<<endl;
    cout<<"imbalance4:"<<imbalance4<<endl;
    cout<<"data_output:" <<data_output<<endl;
    srand(seed);

/* 
int main() 
{



//    string data_input  = "/home/fangjuntao/Brightkite/files/Brightkite/data1,/home/fangjuntao/Brightkite/files/Brightkite/data2";
    //参数
    // string data_input  = "/home/fangjuntao/facebook/data1,/home/fangjuntao/facebook/data2";
    // string data_input = "/home/fangjuntao/Pokec social network/data1,/home/fangjuntao/Pokec social network/data2";
//    string data_input =  "/home/fangjuntao/C++_GraphPartition/C++_GraphPartition/twitch_gamers/data1,/home/fangjuntao/C++_GraphPartition/C++_GraphPartition/twitch_gamers/data2";
    string data_input =  "/home/fangjuntao/kdd2025/LiveJournal/data1,/home/fangjuntao/kdd2025/LiveJournal/data2";
    // string data_input =  "/home/fangjuntao/DBLP/data1,/home/fangjuntao/DBLP/data2";
    
    // unsigned int seed  =  50;

    int  k  =  5;   //k-means 的cluster 数量
    int  p =  10 ;  //partition 的数量
    int    feacuter_size = 3 ;   //注意不同的数据集需要手动设置
    double imbalance1    = 0.10  ; //node
    double imbalance2    = 0.20  ; //feature
    double imbalance3    = 0.70; //innerEdge
    double imbalance4    = 0.70;  //outerEdge
    unsigned int seed    = 50; 
    const int    args_move_sum  = 200000;
    const int    args_swap_sum  = 200000;
    // stringstream(argv[2]) >> imbalance1;
    // stringstream(argv[3]) >> imbalance2;
    // stringstream(argv[4]) >> imbalance3;
    // stringstream(argv[5]) >> imbalance4;
    // stringstream(argv[6]) >> seed;
    // stringstream(argv[7]) >> feacuter_size;
    // stringstream(argv[8]) >> k;，
    double nodeImbalnce       =  numeric_limits<double>::infinity();  //无穷大
    double featureImbalance   =  numeric_limits<double>::infinity();   //无穷大
    double innerEdgeImbance   =  0.0  ; 
    double outerEdgeImbalance =  0.0  ;
    cout<<"seed: "<<seed<<endl;
    cout<<"data_input: "<<data_input<<endl;
    cout<<"imbalance1:"<<imbalance1<<endl;
    cout<<"imbalance2:"<<imbalance2<<endl;
    cout<<"imbalance3:"<<imbalance3<<endl;
    cout<<"imbalance4:"<<imbalance4<<endl;
    srand(seed);
    string data_output = "test.txt";
    // string data_input  =  "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Gowalla/files/Gowalla/data3/part-000.csv";
    // string data_input  = "/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data1/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data2/part-000.csv,/data/home/juntaofang/graphPartitionCopy/GraphPartition_20240223/Brightkite/files/Brightkite/data3/part-000.csv";
    
    //string data_input =  "/data/home/juntaofang/jxqy/data1/part-00000.csv,/data/home/juntaofang/jxqy/data2/part-00000.csv,/data/home/juntaofang/jxqy/data3/part-00000.csv"; 

    
     */
    
    Graph g =  read_graphProcressData(data_input);  //可以优化

    // Do something with the graph
    cout<<" n: "<<g.n<<endl;
    cout<<" m: "<<g.m<<endl; 

    double count_optimal  = g.n/p ;  //最佳点数
    vector<double> vector_optimal(g.feat_sum.size());
    for (size_t i = 0; i < g.feat_sum.size(); ++i) 
    {
        vector_optimal[i] = g.feat_sum[i] /p;
    }

    
    vector<Cluster> result;  //save result  可以优化---> 空间

    vector<int> nodeNumList;
    vector<vector<double>> featuresList;
    vector<int> innerEdgeList ;
    vector<int> outerEdgeList ; 

    double time_start = (double) clock();

    double time_random_start = (double) clock();



/* 

    //随机分配，直到达到要求
    while ( nodeImbalnce > imbalance1 || featureImbalance > imbalance2|| innerEdgeImbance < imbalance3 || outerEdgeImbalance < imbalance4)
    {

        result.clear();
        result.resize(p, Cluster(0,  vector<double>(feacuter_size, 0.0), 0, 0, 0));

        // Iterate over nodes in the graph
        // for (const auto& node_data : g.nodes) 
        for (int i =  0  ; i < g.n ; i++ ) 
        {
            int node = i;

            const  vector<double>& feat = g.nodes[node] ;

            int part_num = rand() % p; // Randomly assign partition
            int j = part_num;
            g.node_id_to_partID[node] = part_num ; // Store partID in the node data

            int new_count = result[j].count + 1;
            for (size_t i = 0; i < result[j].vector.size(); ++i) 
            {
                result[j].vector[i] = (result[j].vector[i] * result[j].count + feat[i]) / new_count;  //存总数
            }
            result[j].count = new_count;
            result[j].point[node]  =  1 ;  //不为0 代表partition包含这个节点 
            result[j].nodeSet.insert(node);  
            int n_neighbors = g.edges[node].size();  //为 0的也会开一个空的unordered_set
            result[j].edgeSum += n_neighbors;  //这里将会得到的是partititon 里面节点的degree和 的两倍

        }


        vector<double> ans  =  computeObject(result,nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );

        nodeImbalnce        =  ans[0];
        double feutureThis = 0.0 ;
        for(int i = 0 ;  i < g.feat_sum.size(); ++i) 
        {
            if(ans[1+i]> feutureThis)  feutureThis =  ans[1+i] ; //存多维特征里面最不均匀的那个特征的imbalance
        }
        featureImbalance    =  feutureThis ; //存多维特征里面最不均匀的那个特征的imbalance
        innerEdgeImbance    =  ans[1+g.feat_sum.size()];
        outerEdgeImbalance  =  ans[2+g.feat_sum.size()]; 



    }


 */


    vector<int>Cluster  =  kmeansGraph(g,k);

    BalancedInitialization(g, Cluster, k , p );

    double time_random_end = (double) clock();

    saveResult(result, g, feacuter_size, p ) ;
    
    vector<double> ans  =  computeObject(result,nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );

    // double time_random_end = (double) clock();
    //k-means先跳过


    //vector<vector<int>> cluster_node_list(1); // 这里由于没有k-means， 所以所有点都是处于同一层   可以优化---> 空间

    size_t i = 0;
    // for (const auto& node_pair : g.nodes) 
    // {
    //     int node = node_pair.first;

    //     cluster_node_list[0].push_back(node);

    //     ++i;
    // }




    //1. 计算gain
    int count_candidate_node;
  //  vector< unordered_map<  string,   vector< pair<int, int>>>> dict_list;
    unordered_map<int,   pair<int,int>> nodeToKey;  //只存gain > 0的点
    unordered_map<int,   int> nodeToGain;   //只存gain > 0的点
    unordered_map<int,  vector<int>> node2neighborsPart;

    double time_gain_start = (double) clock();

    tie(count_candidate_node,  nodeToKey, nodeToGain) = calculateGain( g, node2neighborsPart , p);  //可以优化---> 时间 优化为O(m)
    double time_gain_end = (double) clock();

    cout<<"1125 args_swap_sum: "<<args_swap_sum<<endl;




    //进入move 阶段

    double time_move_start = (double) clock();
    // cout<<" 1435 "<<endl ;
    moveNodesBin(g, nodeToKey, nodeToGain,nodeNumList,featuresList,innerEdgeList,outerEdgeList, node2neighborsPart, p,  args_move_sum, imbalance1, imbalance2, imbalance3,  imbalance4, count_optimal, vector_optimal );
   
    double time_move_end = (double) clock();
    
    double time_end = (double) clock();    
    saveResult(result , g, feacuter_size,p ) ;

    vector<double> ans_move =   computeObject(result, nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );



   
    // 进入swap阶段
    // double time_swap_start = (double) clock();
    // swap_pairs(g , nodeToKey ,cluster_node_list,    dict_list,  
    // node2neighborsPart,  nodeNumList, 
    // featuresList, innerEdgeList,  outerEdgeList, 
    // k,  args_swap_sum, imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  );
 
    
    // saveResult(result ,  g, feacuter_size,k ) ;

    // vector<double> ans_swap =   computeObject(result, nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  k );

    // double time_swap_end = (double) clock();

    
    //把结果写入文件系统

    // string filename = "output.txt";

    write_map_to_file(g,data_output);



    cout << "Algorithm Time: " << (time_end- time_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Balanced Ininital Time: " << (time_random_end- time_random_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Compute Gain  Once Time: " << (time_gain_end- time_gain_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Move Time: " << (time_move_end - time_move_start)/CLOCKS_PER_SEC  << "s" << endl;
    //cout << "Swap Time: " << (time_swap_end - time_swap_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "over!!"<<endl;
    return 0;



}


/* 

    //随机分配，直到达到要求
    while ( nodeImbalnce > imbalance1 || featureImbalance > imbalance2|| innerEdgeImbance < imbalance3 || outerEdgeImbalance < imbalance4)
    {

        result.clear();
        result.resize(p, Cluster(0,  vector<double>(feacuter_size, 0.0), 0, 0, 0));

        // Iterate over nodes in the graph
        // for (const auto& node_data : g.nodes) 
        for (int i =  0  ; i < g.n ; i++ ) 
        {
            int node = i;

            const  vector<double>& feat = g.nodes[node] ;

            int part_num = rand() % p; // Randomly assign partition
            int j = part_num;
            g.node_id_to_partID[node] = part_num ; // Store partID in the node data

            int new_count = result[j].count + 1;
            for (size_t i = 0; i < result[j].vector.size(); ++i) 
            {
                result[j].vector[i] = (result[j].vector[i] * result[j].count + feat[i]) / new_count;  //存总数
            }
            result[j].count = new_count;
            result[j].point[node]  =  1 ;  //不为0 代表partition包含这个节点 
            result[j].nodeSet.insert(node);  
            int n_neighbors = g.edges[node].size();  //为 0的也会开一个空的unordered_set
            result[j].edgeSum += n_neighbors;  //这里将会得到的是partititon 里面节点的degree和 的两倍

        }


        vector<double> ans  =  computeObject(result,nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );

        nodeImbalnce        =  ans[0];
        double feutureThis = 0.0 ;
        for(int i = 0 ;  i < g.feat_sum.size(); ++i) 
        {
            if(ans[1+i]> feutureThis)  feutureThis =  ans[1+i] ; //存多维特征里面最不均匀的那个特征的imbalance
        }
        featureImbalance    =  feutureThis ; //存多维特征里面最不均匀的那个特征的imbalance
        innerEdgeImbance    =  ans[1+g.feat_sum.size()];
        outerEdgeImbalance  =  ans[2+g.feat_sum.size()]; 



    }




 */