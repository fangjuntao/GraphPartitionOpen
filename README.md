# Experimental Environment
The experiments are conducted on a Ubuntu Linux server (Release 22.04.3 LTS) equipped with an Intel Xeon Silver 4210R CPU (10 cores, 20 threads, @ 2.40GHz, 13.8MB L3 Cache) and 512GB memory.
# Summary of Program

# Compile
g++ -o BIHM -std=c++17  -g graph.cpp cluster.cpp evaluation.cpp GainList.cpp BIHM.cpp -O3
# Run
For BIHM, we run the code by inputing the name of the graph file and the budget, e.g., we run the code on Gowalla with budget=100 as follows.

run example：
./BIHM "./Facebook/data1,./Facebook/data2" 0.10 0.20 0.30 0.30 1 5 8 "text.txt"
