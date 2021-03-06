//
// Created by caronkey on 28/12/2021.
//
//

#ifndef SEQGRAPH_MATCHING_H
#define SEQGRAPH_MATCHING_H
#include "Graph.h"
#include "algorithm"
#include <deque>
extern bool VERBOSE;
extern bool BREAK_C;
class matching {
private:
    seqGraph::Graph* graph;
    double** currentMatrix;
    double** originalMatrix;
    std::vector<seqGraph::Vertex*>* originalVertices;
//    std::vector<seqGraph::Junction* >* originalJunctions;
//    seqGraph::Graph* originalGraph;

    int N;

    double* lx;
    double* ly;
    int* matched;

    std::vector<int> cyclePaths;
    void init_labels();
    void update_labels();
    void bfs(int i, double ex[], double ey[], bool visity[],int pre[], double []);
public:
    explicit matching(seqGraph::Graph* graph1);
    ~matching();

    void printM(int i);
    inline int* getMatched() const {
        return matched;
    }
    static bool cmpVertex (int i,int j);
    inline int getN() const {
        return N;
    }
    inline bool isCycle(int i) {
        return std::find(cyclePaths.begin(), cyclePaths.end(), i) != cyclePaths.end();
    }
    inline double** getMatrix() const {
        return this->currentMatrix;
    };
    void hungarian();
    bool dfs(int u, bool visity[], std::vector<int>* pre);
    bool kmDfs(int u, bool visity[],bool visitx[], std::set<int>* pre, double ex[], double ey[], double slack[]);
    void main_steps();

    std::map<int, std::vector<int>*>* resolvePath(std::map<int, std::vector<int>*>* prevPaths);

    int checkConjugateMatch();

//
    std::string idx2StrDir(int idx, const std::string& token= std::string(""));
    std::string idx2Str(int idx);
//    idx to current graph vertex
    seqGraph::Vertex* idx2VertexInCurrentGraph(int idx);
    seqGraph::Vertex* idx2VertexInOriginalGraph(int idx);
    void writeMatchResult(std::ofstream& outS);
    void reconstructMatrix(std::map<int, std::vector<int>*>* paths);
    void resetGraph(seqGraph::Graph* g);
    std::vector<int>* addPrevPath(std::map<int, std::vector<int>*>* prevPath, std::vector<int>* curPath);
    void breakResolvedPaths(std::vector<int>* cur, std::deque<int> & zereBK, std::map<int,std::vector<int>* >* result);
    void breakAndMergeCycle(std::map<int,std::vector<int>*> *result);
    std::vector<int>* breakCycle(std::vector<int> *);
    double* mergePath(std::vector<int>* p1, std::vector<int>* p2, double** matrix, double* result);
};

int conjugateIdx(int idx);
//seqGraph::Graph* reconstructMatrix(double** matrix, std::map<int, std::vector<int>*>*);
std::vector<int>* addPrevPath(std::map<int, std::vector<int>*>* prevPath, std::vector<int>* curPath);
#endif //SEQGRAPH_MATCHING_H
