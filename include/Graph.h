#ifndef __GRAPH_H
#define __GRAPH_H

#include "EventManager.h"
#include <vector>
#include <algorithm>

struct Node;

struct Edge{
    Edge(Node* a, Node* b, EventRef eRef):
        ref(eRef), first(a), second(b)
    {}
    EventRef ref;
    Node *first, *second;
};

struct Node{
    std::vector<Edge*> edges;
};

class Graph{
public:
    Graph(size_t nEntries):
        N_(nEntries)
    {
        nodes_ = new Node[N_];
    }
    ~Graph(void){
        for(size_t i = 0; i < N_; ++i){
            clear(i);
        }
        delete nodes_;
    }
    void addEdge(size_t a, size_t b, EventRef eRef){
        Edge* edge = new Edge(&nodes_[a], &nodes_[b], eRef);
        nodes_[a].edges.push_back(edge);
        nodes_[b].edges.push_back(edge);
    }
    std::vector<EventRef> getAssociations(size_t a){
        std::vector<EventRef> result;
        for(auto edge: nodes_[a].edges){
            result.push_back(edge->ref);
        }
        return result;
    }
    void clear(size_t a){
        for(auto edge: nodes_[a].edges){ //Loop through all the edges of the node to be cleared
            Node* otherNode = (&nodes_[a] != edge->first)? edge->first: edge->second; //get the node on the other side of edge
            otherNode->edges.erase(  //Erase edge in other node's edge list
                std::find_if(otherNode->edges.begin(), otherNode->edges.end(),
                    [&edge](const Edge* nedge){
                        return nedge == edge;
                    }
                )
            );
            delete edge;
        }
        nodes_[a].edges.clear();
    }
private:
    size_t N_;
    Node *nodes_;
};

#endif
