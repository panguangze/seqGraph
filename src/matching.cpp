//
// Created by caronkey on 28/12/2021.
//

#include "../include/matching.h"
#include <cstring>
#include <algorithm>
#include <queue>
#include <iostream>
const double ZERO = 0.00000001;

void matching::printM(int i){
    std::cout<<i<<"|";
    for (int j = 0; j <  N; j++){
        std::cout<<this->getMatrix()[i][j]<<"\t";
    }
}
const double INF = 1e18;
matching::matching(seqGraph::Graph* graph1) {
    this->graph = graph1;
    N = 2 * graph1->getVCount();
    this->matched = new int[N + 1];
    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
}

bool matching::kmDfs(int u, bool visity[],bool visitx[], std::vector<int>* pre, double ex[], double ey[], double slack[]) {
    auto matrix = getMatrix();
    visity[conjugateIdx(u)] = true;
    visitx[u] = true;
    for (int i = 1; i < N+1; i++) {
        if(i == u) continue;
//        auto t1 = ex[u];
//        auto t2 = ey[i];
//        auto t3 = matrix[u][i];
//        auto t4 = visity[i];
//        auto t5 = matched[i];
//        auto t6 = slack[i];
//        auto p = std::find(pre->begin(), pre->end(), i) != pre->end();
        if (visity[i]) continue;
//        if ( u == 19 && i == 1615) {
//            std::cout<<matched[i];
//            printM(u);
//            std::cout<<std::endl;
//             int km = 0;
//            std::cout<<std::endl;
//        }
        if(std::abs(ex[u] + ey[i] - matrix[u][i]) <= ZERO) {
            visity[i] = true;
            if((matched[i] == -1) || (kmDfs(matched[i], visity, visitx, pre, ex, ey, slack))) {
                pre->push_back(u);
                pre->push_back(conjugateIdx(i));
                matched[i] = u;
                matched[conjugateIdx(u)] = conjugateIdx(i);
                return true;
            }
        } else if(slack[i] > ex[u] + ey[i] - matrix[u][i]) {
//            auto tt1 = ex[u];
//            auto tt2 = ey[i];
//            auto tt3 = matrix[u][i];
            slack[i] = ex[u] + ey[i]  - matrix[u][i];
//            auto tts = slack[i];
//            auto ttttt = 99;
        }
    }
    return false;
}


bool matching::dfs(int u, bool visity[], std::vector<int>* pre) {
    auto matrix = getMatrix();
    visity[conjugateIdx(u)] = true;
    for (int i = 1; i < N+1; i++) {
//        auto p = std::find(pre->begin(), pre->end(), i) != pre->end();
        if(matrix[u][i] != 0 && !visity[i]) {
            visity[i] = true;
            if((matched[i] == 0) || (dfs(matched[i], visity, pre))) {
                pre->push_back(u);
                pre->push_back(conjugateIdx(i));
                matched[i] = u;
                matched[conjugateIdx(u)] = conjugateIdx(i);
                return true;
            }
        }
    }
    return false;
}

void matching::bfs(int u, double ex[], double ey[], bool visity[], int pre[], double slack[]) {
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
        visity[conjugateIdx(x)] = true;
        for (int i = 1; i < N + 1; i++) {
            auto cI = conjugateIdx(i);
            if(!visity[i]) {
                if(slack[i] > ex[x] + ey[i] - matrix[x][i]){
                    slack[i] = ex[x] + ey[i] - matrix[x][i];
                    pre[i] = y;
//                    pre[conjugateIdx(y)] = cI;
                }
                if(slack[i] < delta) {
                    delta = slack[i], yy = i;
                }
            }
        }
        for(int i = 0; i < N + 1; i++){
            if(visity[i]){
                ex[matched[i]] -= delta;
                ey[i] += delta;
//                ex[conjugateIdx(i)] -= delta;
//                ey[conjugateIdx(matched[i])] += delta;
            }
            else slack[i] -= delta;
        }
        y = yy;
    }
    while(y){
        matched[y] = matched[pre[y]];
        auto t1 = conjugateIdx(matched[pre[y]]);
        auto t2 = conjugateIdx(y);
        matched[t1] = t2;
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

void matching::hungarian() {
    auto matrix = getMatrix();
    bool visity[N + 1];
    bool visitx[N + 1];
    double slack[N + 1], ex[N + 1], ey[N + 1];

    memset(ey,0,sizeof ey);
    memset(visity,0,sizeof visity);

    for( int i = 1 ; i <= N ; i++ ){
        ex[i] = -INF;
        for( int j = 1 ; j <= N ; j++ ){
            if( ex[i] < matrix[i][j] ) ex[i] = matrix[i][j];
        }
    }
    auto* pre = new std::vector<int>();
    for (int i = 1; i < N + 1; i++) {
        for( int l = 1 ; l <=N ; l++ ) slack[l] = 1000;
        if (std::find(pre->begin(), pre->end(), i) != pre->end()) {
            continue;
        }
        while (true) {
            memset(visity, 0, sizeof visity);
            memset(visitx, 0, sizeof visitx);
            auto r = kmDfs(i, visity,visitx, pre,ex, ey, slack);
            if (r) {
                break;
            } else {
                double delta = INF;
                for( int j = 1 ; j <= N ; j++ ){
                    if( !visity[j] && delta > slack[j] ){
                        delta = slack[j];
                    }
                }

                for( int k = 1 ; k <= N; k++ ){
                    if( visitx[k] ) ex[k] -= delta;
                    if( visity[k] ) ey[k] += delta;
                    else slack[k] -= delta;
                }
            }
        }
        if (VERBOSE) {
            std::cout<<i<<"\t|";
            for(int k = 0; k < N + 1; k++) {
                std::cout<<matched[k]<<"\t";
            }
            std::cout<<std::endl;
        }
    }
}
void matching::main_steps() {
    bool visity[N + 1];
    int pre[N + 1];
    double slack[N + 1], ex[N + 1], ey[N + 1];

    memset(ex,0,sizeof ex);
    memset(ey,0,sizeof ey);
    memset(visity,0,sizeof visity);
    memset(pre, 0, sizeof pre);
    std::vector<int> skipped;
    for(int i = 1; i < N + 1; ++ i){
//        if (i % 2 != 1) continue;
        for(int j = 1; j < N + 1; ++ j)
            ey[j] = false;
        if(std::find(skipped.begin(), skipped.end(), i) != skipped.end()) continue;
        bfs(i, ex, ey, visity, pre, slack);
        auto cI = conjugateIdx(i);
        skipped.push_back(matched[cI]);
        for(int k = 0; k < N + 1; k++) {
            std::cout<<matched[k]<<"\t";
        }
        std::cout<<std::endl;
    }
}

std::map<int, std::vector<std::string>*>* matching::resolvePath() {
    auto matrix = graph->getConjugateMatrix();
    bool visited[N+1];
    memset(visited, 0, sizeof visited);
    auto* resolvePath = new std::map<int, std::vector<std::string> *>();
    auto checkedC = checkConjugateMatch();
    if(checkedC!= 0)
        std::cout<<"Conjugate not checked\n";
    else
        std::cout<<"check conjugate done\n";
    std::cout<<"final paths\n";
    for(int i = 1; i < N + 1; i++) {
        if (visited[i]) continue;
        if (matrix[matched[i]][i] == 0) continue;
        auto* currentPath = new std::vector<std::string>();
        int now = i;
        int vIdx = (now + 1) / 2;
        char dir = now % 2 == 0 ? '-':'+';
        currentPath->push_back((*this->graph->getVertices())[vIdx - 1]->getId() + dir);
        if ((*this->graph->getVertices())[vIdx - 1]->getId() == "52_0")
            int tm = 99;

//        std::cout<< (*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
        visited[now] = true;
        visited[conjugateIdx(now)] = true;
        bool currentInsert = true;
        while (matrix[matched[now]][now] != double(0)) {
            if (resolvePath->find(matched[now]) != resolvePath->end()) {
                auto prevPath = (*resolvePath)[matched[now]];
//                +"_"+ std::to_string(now)
                for (auto it = currentPath->rbegin(); it != currentPath->rend(); it++){
                    prevPath->insert(prevPath->begin(), *it);
                }
                resolvePath->erase(matched[now]);
                resolvePath->emplace(i, prevPath);
                currentInsert = false;
                break;
            }
            if(visited[matched[now]]){
                if (matched[now] == i) currentPath->push_back("c");
                break;
            }
            visited[matched[now]] = true;
            visited[conjugateIdx(matched[now])] = true;
            vIdx = (matched[now] + 1) / 2;
            dir = matched[now] % 2 == 0 ? '-':'+';
//            std::cout<<(*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
            now = matched[now];
            currentPath->push_back((*this->graph->getVertices())[vIdx - 1]->getId() + dir);
        }
        if(currentInsert)
            (*resolvePath)[i] = currentPath;
    }
    return resolvePath;
}

int conjugateIdx(int idx) {
    int cI;
    if(idx == 0) return idx;
    cI = idx%2 == 0? idx-1:idx+1;
    return cI;
}

int matching::checkConjugateMatch() {
    int r = 0;
    auto matrix = this->getMatrix();
    for (int i = 1; i < N+1; i++) {
        auto left = matched[i];
        if (matrix[left][i] == 0) continue;
        auto conjugateI = conjugateIdx(i);
        auto conjugateLeft = conjugateIdx(left);
        if (matched[conjugateLeft] != conjugateI) {
            r++;
            std::cout<<i<<'\t'<<matched[i]<<'\t'<<conjugateLeft<<'\t'<<conjugateI<<'\t'<<matched[conjugateLeft]<<'\n';
        }
    }
    std::cout<<std::endl;
    return r;
}