//
// Created by caronkey on 28/12/2021.
//
//

#ifndef SEQGRAPH_MATCHING_H
#define SEQGRAPH_MATCHING_H
#include "Graph.h"
class matching {
private:
    seqGraph::Graph* graph;

    int N;

    double* lx;
    double* ly;
    int* matched;
    void init_labels();
    void update_labels();
    void bfs(int i, double* ex, double* ey, bool* visity,int* pre, double* slack);
public:
    explicit matching(seqGraph::Graph* graph1);
    inline int* getMatched() const {
        return matched;
    }
    inline int getN() const {
        return N;
    }
    inline double** getMatrix() const {
        return this->graph->getConjugateMatrix();
    };
    ~matching() {
        free(this->graph);
    }
    void hungarian();
    void main_steps();

    void resolvePath();
};

int conjugateIdx(int idx);

#endif //SEQGRAPH_MATCHING_H
