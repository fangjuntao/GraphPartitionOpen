 g++ -o BIHM -std=c++17  -g graph.cpp cluster.cpp evaluation.cpp GainList.cpp BIHM.cpp -O3


run example：
./BIHM "./Facebook/data1,./Facebook/data2" 0.10 0.20 0.30 0.30 1 5 8 "text.txt"
