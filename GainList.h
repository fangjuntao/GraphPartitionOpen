#ifndef GAINLIST_H
#define GAINLIST_H

#include <vector>
#include <map>
#include <list>
#include <unordered_map>

struct Node {
    int id;
    int gain;
    
    Node() {}

    Node(int id, int gain) : id(id), gain(gain) {}
};

class GainList {
public:
    GainList(int size);
    void insert(Node node);
    void remove(int nodeId);
    void update(int nodeId, int newGain);
    void print();
    int size() const;
    const std::list<Node>& operator[](int index) const;
    void clear();
    void resize(int newSize);

private:
    std::vector<std::list<Node>> gainLists;
    std::unordered_map<int, std::list<Node>::iterator> nodePositions;
};

#endif // GAINLIST_H