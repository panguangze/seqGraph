//
// Created by caronkey on 28/12/2021.
//

#include "../include/matching.h"
#include <cstring>
#include <algorithm>
#include <queue>
#include <iostream>

const double INF = 1e18;
matching::matching(seqGraph::Graph* graph1) {
    this->graph = graph1;
    N = 2 * graph1->getVCount();
    this->matched = new int[N + 1];
}

void matching::bfs(int u, double * ex, double* ey, bool* visity, int* pre, double* slack) {
    auto matrix = this->graph->getConjugateMatrix();
    int x,cY,y=0, yy=0;
    double delta;

    for (int i = 1; i < N + 1; i++) slack[i] = INF;

//    当前的match以及共轭的match
    matched[y] = u;
//    matched[conjugateIdx(u)] = y;
    while (matched[y]) {
        x = matched[y], delta = INF;
//        cY = matched[conjugateIdx(u)];
//      当前点y和x共轭点同时被visit
        visity[y] = true;
//        visity[conjugateIdx(x)] = true;
        for (int i = 1; i < N + 1; i++) {
            if(!visity[i]) { // visty and its conjugate will always sync.
                if(slack[i] > ex[x] + ey[i] - matrix[x][i]){
                    slack[i] = ex[x] + ey[i] - matrix[x][i];
                    pre[i] = y;
                }
                if(slack[i] < delta) {
                    delta = slack[i], yy = i;
                }
            }
        }
        for(int i = 0; i < N + 1; i++){
            if(visity[i]) ex[matched[i]] -= delta,ey[i] += delta;
            else slack[i] -= delta;
        }
        y = yy;
    }
    while(y){
        matched[y] = matched[pre[y]];
//        auto t1 = conjugateIdx(matched[pre[y]]);
//        auto t2 = conjugateIdx(y);
//        matched[t1] = conjugateIdx(y);
        y = pre[y];
    }
}


//void matching::bfs(int u, double * ex, double* ey, bool* visity, int* pre, double* slack) {
//    auto matrix = this->graph->getConjugateMatrix();
//    int x,cY,y=0, yy=0;
//    double delta;
//
//    for (int i = 1; i < N + 1; i++) slack[i] = INF;
//
////    当前的match以及共轭的match
//    matched[y] = u;
////    matched[conjugateIdx(u)] = y;
//    while (matched[y]) {
//        x = matched[y], delta = INF;
////        cY = matched[conjugateIdx(u)];
////      当前点y和x共轭点同时被visit
//        visity[y] = true;
////        visity[conjugateIdx(x)] = true;
//        for (int i = 1; i < N + 1; i++) {
//            if(!visity[i]) { // visty and its conjugate will always sync.
//                if(slack[i] > ex[x] + ey[i] - matrix[x][i]){
//                    slack[i] = ex[x] + ey[i] - matrix[x][i];
//                    pre[i] = y;
//                }
//                if(slack[i] < delta) {
//                    delta = slack[i], yy = i;
//                }
//            }
//        }
//        for(int i = 0; i < N + 1; i++){
//            if(visity[i]) ex[matched[i]] -= delta,ey[i] += delta;
//            else slack[i] -= delta;
//        }
//        y = yy;
//    }
//    while(y){
//        matched[y] = matched[pre[y]];
////        auto t1 = conjugateIdx(matched[pre[y]]);
////        auto t2 = conjugateIdx(y);
////        matched[t1] = conjugateIdx(y);
//        y = pre[y];
//    }
//}
void matching::main_steps() {
    bool visity[N + 1];
    int pre[N + 1];
    double slack[N + 1], ex[N + 1], ey[N + 1];

    memset(ex,0,sizeof ex);
    memset(ey,0,sizeof ey);
    memset(visity,0,sizeof visity);
    memset(pre, 0, sizeof pre);

    for(int i = 1; i < N + 1; ++ i){
        for(int j = 1; j < N + 1; ++ j)
            ey[j] = false;
        bfs(i, ex, ey, visity, pre, slack);
        for(int k = 0; k < N + 1; k++) {
            std::cout<<matched[k]<<"\t";
        }
        std::cout<<std::endl;
    }
}

int conjugateIdx(int idx) {
    int cI;
    if(idx == 0) return idx;
    cI = idx%2 == 0? idx-1:idx+1;
    return cI;
}

void matching::resolvePath() {
    auto matrix = graph->getConjugateMatrix();
    bool visited[N+1];
    memset(visited, 0, sizeof visited);
    for(int i = 1; i < N + 1; i++) {
        if (visited[i]) continue;
        if (matrix[matched[i]][i] == 0) continue;
        int now = i;
        int vIdx = (now + 1) / 2;
        char dir = now % 2 == 0 ? '-':'+';
        std::cout<< (*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
        visited[now] = true;
        visited[conjugateIdx(now)] = true;
        while (matrix[matched[now]][now] != double(0)) {
            if(visited[matched[now]]){
                std::cout<<'c';
                break;
            }
            visited[matched[now]] = true;
            visited[conjugateIdx(matched[now])] = true;
            vIdx = (matched[now] + 1) / 2;
            dir = matched[now] % 2 == 0 ? '-':'+';
            std::cout<<(*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
            now = matched[now];
        }
        std::cout<<'\n';
    }
}