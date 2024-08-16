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



using FourTuple = std::tuple<int, int, int, int>;
using FiveTuple = std::tuple<bool, int, int, int, int>;


// Function to normalize the features of nodes
vector<vector<double>> normalize_features(const std::vector<std::vector<double>>& nodes_orign) 
{

    vector<vector<double>> nodes  = nodes_orign ;
    int nobs = nodes.size();  
    int ndim = nodes[0].size();  

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




    
    vector<vector<double>> nodes  =  normalize_features(g.nodes);

    int nobs = nodes.size();  
    int ndim = nodes[0].size();  
    int ncenters = k;  


    std::vector<double> matrix;
    for (const auto& node : nodes) 
    {
        matrix.insert(matrix.end(), node.begin(), node.end());
    }

    kmeans::SimpleMatrix<double, int, int> kmat(ndim, nobs, matrix.data());


    auto res = kmeans::compute(
        kmat,
        kmeans::InitializeKmeanspp<kmeans::SimpleMatrix<double, int, int>, int, double>(),  
        kmeans::RefineLloyd<kmeans::SimpleMatrix<double, int, int>, int, double>(), 
        ncenters
    );

    vector<int> ans  = res.clusters; 


 
    return ans;

}


void BalancedInitialization(Graph &g, const vector<int> & kmeansResult, int  k , int p ) //均匀初始化
{
    vector<int> clusterCount(k,0);
    unordered_map<int,bool> vertexState;  
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


    std::vector<std::pair<int, int>> nodeNeighborCount;
    for (int i = 0; i < g.edges.size(); ++i) {
        nodeNeighborCount.push_back({i, g.edges[i].size()});
    }


    std::sort(nodeNeighborCount.begin(), nodeNeighborCount.end(),
            [](const std::pair<int, int>& a, const std::pair<int, int>& b) 
            {
                return a.second > b.second;
            });

 
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

    

        auto minElementIt  = std::min_element(score.begin(), score.end());


        int selectedPartID = std::distance(score.begin(), minElementIt);

        g.node_id_to_partID[node]   =  selectedPartID ;




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



void moveNodesBin(Graph &g,  unordered_map<int,   pair<int,int> >& nodeToKey , unordered_map<int,  int>& nodeToGain, vector<int>& nodeNumList, vector<vector<double>>& featuresList, vector<int>& innerEdgeList, vector<int>& outerEdgeList, unordered_map<int,  vector<int>> &node2neighborsPart, int k,  int args_swap_sum, double imbalance1, double imbalance2, double imbalance3, double imbalance4, int count_optimal, const vector<double>& vector_optimal ) //省空间版本，加速
{
  
    unordered_map<int, FourTuple> node_FourTuple ;  // （node, value1,value2,value3,value4)


   // unordered_map<int,  int> nodeToGain ;   //nodeID -gain
    unordered_map<int, bool> nodeToConstraints ;  // nodeID--bool  

    vector<pair<int,int>> insertVector;  //index, value
    vector<pair<int,int>> deleteVector ;  //index,value




    GainList gainList(10000);  
    
    int swap_real_num = 0; 

    // vector<int> max_gain_index ;
    int max_gain_sum = 0 ; 

    int max_gain               = 0 ;
    int count_ban_innerEdge    = 0 ;
    int count_ban_outerEdge    = 0 ;
    int count_ban_nodeNum      = 0 ;
    int count_ban_features     = 0 ;
    int count_candidate_node   = 0 ;

    int count_modifiNodes =  0;

    int index_1_out_value      = -1 ; 
    int index_2_out_value      = -1 ;
    int index_1_inner_value    = -1 ;
    int index_2_inner_value    = -1 ; 

    for(int i = 0 ; i< g.n ; i++)   
    {
        int  node        = i;
        int  gain        =  0 ; 
        if( nodeToGain.find(node) != nodeToGain.end() )
        {
            gain = nodeToGain.at(node);  
        }
        else
        {
            continue;
        }





        int index_1 = nodeToKey.at(node).first ;
        
        int index_2 = nodeToKey.at(node).second;
       // getIndex(key,index_1,index_2) ;   // Convert string to int  pair<int, int>
        
        
        
        // auto node_gain   = value[i];

        
        bool  constraints ; 
        tie(constraints, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
            IfSatisfiedContrainst(g, node, index_1,  index_2, node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

        node_FourTuple[node]  =  make_tuple(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value);  
        //bool constraints  =   IfSatisfiedContrainst(node, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;
        

       // nodeToGain[node]   =   gain  ;

        if( !constraints )  
        {
            // nodeToGain[node]  =  -1  ;  
            continue;
        }
        else 
        {
            nodeToConstraints[node]  =  constraints ;  

            if (gain > 0 )  
            {


                gainList.insert(Node(node,gain));
                max_gain   = max(max_gain, gain);
                ++count_candidate_node ;
            }
        
        }


    
    
    
    }

    // cout<< "max_gain: "<< max_gain<<endl;
    
    

    int findMaxNum  =  0 ;
    while(1)       
    {
        int swap_real_numCopy  =  swap_real_num;

        for(int i = max_gain; i > 0; --i )    
        {
            deleteVector.clear();
            insertVector.clear();
         
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

                    int index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value;
                    tie(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = node_FourTuple.at(node_1.id);



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


                    int returnValue1  = nodeToKey.erase(node_1.id) ;
                    int returnValue2 =  nodeToGain.erase(node_1.id);
                    if(returnValue1 == 0  || returnValue2 == 0 )
                    {
                        cout<<" wrong in main line 245***************"<<endl;
                    }



                    deleteVector.push_back(make_pair(gain, node_1.id));
                    nodeToGain.erase(node_1.id);
                    nodeToKey.erase(node_1.id);
                    nodeToConstraints.erase(node_1.id);


                    for(const auto & neighbor: g.edges.at(node_1.id))
                    {

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
                    



                        if ( nodeToGain.find(neighbor) != nodeToGain.end() && nodeToConstraints.find(neighbor) !=nodeToConstraints.end() )  //之前在vector里面 
                        {
                            




                            if(nodeToGain[neighbor] != gain)
                            {


   
                                if(gain == 0 ||constraints== false ) 
                                {


                                    deleteVector.push_back(make_pair(nodeToGain[neighbor],neighbor));

                                    
                                    --count_candidate_node ;

                                
                                }
                                else 
                                {


                                    deleteVector.push_back(make_pair(nodeToGain[neighbor],neighbor));

                                    insertVector.push_back(make_pair(gain, neighbor));
                                    max_gain  =  max(gain, max_gain);
                            
                                
                                }

                            }
                            else  
                            {
            

                                if( constraints != nodeToConstraints[neighbor])
                                {


                                    deleteVector.push_back(make_pair(nodeToGain[neighbor], neighbor) );
                                    --count_candidate_node ;

                                }


                            }

                        }
                        else  
                        {


                            if( gain > 0 && constraints)
                            {


                                
                                //constraintsNodeSet[gain].insert(neighbor);
                                insertVector.push_back(make_pair(gain,neighbor ));
                                max_gain  =  max(gain, max_gain);
                                count_candidate_node ++;
                            }



                        

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

                  

                    max_gain_sum += max_gain;


                    if (swap_real_num % 100 == 0)
                    {
                        
                        cout<<"swap_real_num: "<<swap_real_num<<endl;
                        

                        auto current_time =    chrono::system_clock::now();
                        

                        time_t time_t_current_time =    chrono::system_clock::to_time_t(current_time);


                        string time_str =    ctime(&time_t_current_time);
                        


                        cout <<"max_gain:"<<max_gain<<endl;

                        cout<<" max_gain_sum :"<<max_gain_sum<<endl;
                        cout<<"count_ban_innerEdge:"<<count_ban_innerEdge<<endl;
                        cout<<"count_ban_outerEdge:"<<count_ban_outerEdge<<endl;
                        cout<<"count_ban_feauture:"<< count_ban_features<<endl;
                        cout<<"count_ban_nodeNum: "<<count_ban_nodeNum<<endl;
                        cout<<"count_candidate_node: "<<count_candidate_node<<endl;
                        cout<<" count_modifiNodes: "<<count_modifiNodes<<endl;
                        // cout<<" 100 轮总共查找最大值的次数："<<findMaxNum<<endl;
                        cout << "Current time: " << time_str << endl;
                        findMaxNum  = 0 ;
                        count_modifiNodes = 0 ;


                    }


                    i =  0 ;
                }



                break;
            }

            count_modifiNodes += deleteVector.size();
            for(auto &pair: deleteVector)
            {

                int index =  pair.first;
                int node  =  pair.second;

                gainList.remove(node);


            }
            count_modifiNodes += insertVector.size();
            for(auto& pair:insertVector)
            {

                int index =  pair.first;
                int node  =  pair.second;


                gainList.insert(Node(node,index));
                   


            }
        
        }

        if(swap_real_numCopy == swap_real_num)
        {
            // constraintsNodeSet.clear();
            // constraintsNodeSet.resize(g.n);
            gainList.clear();
            gainList.resize(g.n);
            max_gain =  0;
            count_candidate_node  = 0 ;

                    
            for(int i  =0 ; i< g.n ; i++)   //O(n)    
            {
                int  node        =  i ;
                int  gain        =  0 ; 
                if( nodeToGain.find(node) != nodeToGain.end() )
                {
                    gain = nodeToGain.at(node);  
                }
                else
                {
                    continue;
                }



                int index_1 = nodeToKey.at(node).first ;
                
                int index_2 = nodeToKey.at(node).second;
                // getIndex(key,index_1,index_2) ;   // Convert string to int  pair<int, int>
                
                
                
                // auto node_gain   = value[i];

                
                bool  constraints ; 
                tie(constraints, index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value) = 
                    IfSatisfiedContrainst(g, node, index_1,  index_2,node2neighborsPart, nodeNumList, featuresList, innerEdgeList, outerEdgeList,imbalance1, imbalance2, imbalance3, imbalance4, count_optimal,vector_optimal  ) ;

                node_FourTuple[node]  =  make_tuple(index_1_out_value , index_2_out_value, index_1_inner_value, index_2_inner_value);  

                if( !constraints )
                {
                    
                    continue;
                }
                else
                {
                    nodeToConstraints[node]  =  constraints ;   

                    if (gain > 0 )  //gain 一定大于0
                    {


                        gainList.insert(Node(node,gain));
                        max_gain   = max(max_gain, gain);
                        ++count_candidate_node ;
                    }
                
                }


            
            
            
            }


            
            if(max_gain == 0)
            {

                //cout<<"无法继续move"<<endl;
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



int main(int argc, char* argv[]) 
{

    if(argc < 10) 
    {
        cout << "please input: data_input, imbalance1, imbalance2, imbalance3, imbalance4, feacuter_size, k, p, output_file" << endl;
        return 1;
    }
    string data_input  = argv[1];
    string data_output = argv[9];
    


    int  k, p ;  
    int    feacuter_size  ;  
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
    // stringstream(argv[6]) >> seed;
    stringstream(argv[6]) >> feacuter_size;
    stringstream(argv[7]) >> k;
    stringstream(argv[8]) >> p;
    double nodeImbalnce       =  numeric_limits<double>::infinity();  
    double featureImbalance   =  numeric_limits<double>::infinity();   
    double innerEdgeImbance   =  0.0  ; 
    double outerEdgeImbalance =  0.0  ;
    // cout<<"seed: "<<seed<<endl;
    seed = 50;
    cout<<"data_input: "<<data_input<<endl;
    cout<<"imbalance1:"<<imbalance1<<endl;
    cout<<"imbalance2:"<<imbalance2<<endl;
    cout<<"imbalance3:"<<imbalance3<<endl;
    cout<<"imbalance4:"<<imbalance4<<endl;
    cout<<"data_output:" <<data_output<<endl;
    imbalance3 =  1-imbalance3;
    imbalance4 =  1-imbalance4;
    srand(seed);
    Graph g =  read_graphProcressData(data_input);  
    // Do something with the graph
    cout<<" n: "<<g.n<<endl;
    cout<<" m: "<<g.m<<endl; 

    double count_optimal  = g.n/p ;  
    vector<double> vector_optimal(g.feat_sum.size());
    for (size_t i = 0; i < g.feat_sum.size(); ++i) 
    {
        vector_optimal[i] = g.feat_sum[i] /p;
    }

    
    vector<Cluster> result;  //save result  

    vector<int> nodeNumList;
    vector<vector<double>> featuresList;
    vector<int> innerEdgeList ;
    vector<int> outerEdgeList ; 

    double time_start = (double) clock();

    double time_random_start = (double) clock();


    vector<int>Cluster  =  kmeansGraph(g,k);

    BalancedInitialization(g, Cluster, k , p );

    double time_random_end = (double) clock();

    saveResult(result, g, feacuter_size, p ) ;
    
    vector<double> ans  =  computeObject(result,nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );





    size_t i = 0;






    int count_candidate_node;
  //  vector< unordered_map<  string,   vector< pair<int, int>>>> dict_list;
    unordered_map<int,   pair<int,int>> nodeToKey;  
    unordered_map<int,   int> nodeToGain;  
    unordered_map<int,  vector<int>> node2neighborsPart;

    double time_gain_start = (double) clock();

    tie(count_candidate_node,  nodeToKey, nodeToGain) = calculateGain( g, node2neighborsPart , p);  //可以优化---> 时间 优化为O(m)
    double time_gain_end = (double) clock();

    cout<<"1125 args_swap_sum: "<<args_swap_sum<<endl;






    double time_move_start = (double) clock();

    moveNodesBin(g, nodeToKey, nodeToGain,nodeNumList,featuresList,innerEdgeList,outerEdgeList, node2neighborsPart, p,  args_move_sum, imbalance1, imbalance2, imbalance3,  imbalance4, count_optimal, vector_optimal );
   
    double time_move_end = (double) clock();
    
    double time_end = (double) clock();    
    saveResult(result , g, feacuter_size,p ) ;

    vector<double> ans_move =   computeObject(result, nodeNumList,featuresList, innerEdgeList,outerEdgeList,  g, count_optimal, vector_optimal,  p );




    write_map_to_file(g,data_output);



    cout << "Algorithm Time: " << (time_end- time_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Balanced Ininital Time: " << (time_random_end- time_random_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Compute Gain  Once Time: " << (time_gain_end- time_gain_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "Move Time: " << (time_move_end - time_move_start)/CLOCKS_PER_SEC  << "s" << endl;
    //cout << "Swap Time: " << (time_swap_end - time_swap_start)/CLOCKS_PER_SEC  << "s" << endl;
    cout << "over!!"<<endl;
    return 0;



}

