# Experimental Environment
The experiments are conducted on an Ubuntu Linux server (Release 22.04.3 LTS) equipped with:

* CPU: Intel Xeon Silver 4210R (10 cores, 20 threads, @ 2.40GHz, 13.8MB L3 Cache)
* Memory: 512GB
# Summary of Program
This project consists of several C++ files:

* graph.cpp: Contains the data structure for storing the attribute graph.
* cluster.cpp: Contains the data structure for storing information about every partition.
* evaluation.cpp:  Calculates various metrics for evaluating the graph partitioning results.
* GainList.cpp: Implements the $GList$ data structure as described in the paper.
* BIHM.cpp: Implements the BIHM algorithm.

# Compile
To compile the program, use the following command:

```
g++ -o BIHM -std=c++17  -g graph.cpp cluster.cpp evaluation.cpp GainList.cpp BIHM.cpp -O3
```

# Execution Instructions
To run the BIHM program, use the following command format:
```
./BIHM <data_input> <imbalance1> <imbalance2> <imbalance3> <imbalance4> <feature_size> <k> <p> <out_putfile>
```

## Input Parameter
<data_input>: The path to the input graph files. The graph dataset is divided into two parts:
1. Path to the graph structure information (e.g., ./Facebook/data1).
2. Path to the graph attribute information (e.g., ./Facebook/data2).

`<imbalance1>`: the balance parameters of vertex

`<imbalance2>`: the balance parameters of attribute

`<imbalance3>`: the balance parameters of edge

`<imbalance4>`: the balance parameters of cut

`<feature_size>`: number of user attributes

`<k>`: the $k$ of $k$-means

`<p>`: the number of partitions 
`<out_putfile>`: the path to the output file where the graph partitioning results will be saved.

# Example Usage
Here is an example of how to run the program:
```
./BIHM "./Facebook/data1,./Facebook/data2" 0.10 0.20 0.30 0.30 1 5 8 "results.txt"
```
This command will execute the BIHM program using the specified input datasets (Facebook) and parameters to partition the Facebook Graph into 8 paritions, and the partitioning results will be saved in results.txt.
