//
// Created by caronkey on 28/12/2021.
//
//

#ifndef SEQGRAPH_MATCHING_H
#define SEQGRAPH_MATCHING_H
#include "Graph.h"
extern int VERBOSE;
class matching {
private:
    seqGraph::Graph* graph;

    int N;

    double* lx;
    double* ly;
    int* matched;
    void init_labels();
    void update_labels();
    void bfs(int i, double ex[], double ey[], bool visity[],int pre[], double []);
public:
    explicit matching(seqGraph::Graph* graph1);
    void printM(int i);
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
    bool dfs(int u, bool visity[], std::vector<int>* pre);
    bool kmDfs(int u, bool visity[],bool visitx[], std::vector<int>* pre, double ex[], double ey[], double slack[]);
    void main_steps();

    std::map<int, std::vector<std::string>*>* resolvePath();

    int checkConjugateMatch();

};

int conjugateIdx(int idx);

#endif //SEQGRAPH_MATCHING_H
