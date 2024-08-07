#include "GainList.h"
#include <iostream>

using namespace std;

GainList::GainList(int size) {
    gainLists.resize(size);
}

void GainList::insert(Node node) {
    int gain = node.gain;
    if (gain >= (int)gainLists.size()) {
        gainLists.resize(gain + 1);
    }
    gainLists[gain].push_front(node);
    nodePositions[node.id] = gainLists[gain].begin();
}

void GainList::remove(int nodeId) {
    auto it = nodePositions.find(nodeId);
    if (it != nodePositions.end()) {
        int gain = it->second->gain;
        gainLists[gain].erase(it->second);
        nodePositions.erase(it);
    }
    else{
       std::cout<<"Gainlist.cpp line 27: erro"<<std::endl;
    }
}

void GainList::update(int nodeId, int newGain) {
    remove(nodeId);
    Node newNode;
    newNode.id = nodeId;
    newNode.gain = newGain;
    insert(newNode);
}

void GainList::print() {
    for (int i = 0; i < gainLists.size(); ++i) {
        cout << "Gain " << i << ": ";
        for (const auto& node : gainLists[i]) {
            cout << node.id << " ";
        }
        cout << endl;
    }
}

int GainList::size() const {
    return gainLists.size();
}

const std::list<Node>& GainList::operator[](int index) const {
    return gainLists[index];
}

void GainList::clear() {
    for (auto& gainList : gainLists) {
        gainList.clear();
    }
    nodePositions.clear();
}

void GainList::resize(int newSize) {
    gainLists.resize(newSize);
}