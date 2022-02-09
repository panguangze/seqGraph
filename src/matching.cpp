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
    std::fill_n(this->matched,N+1, -1);
//    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
    currentMatrix = this->graph->getConjugateMatrix();
//    this->originalGraph = graph1;
    this->originalMatrix = new double*[N+1];
    for(int i = 0; i < N+1; ++i)
        originalMatrix[i] = new double [N+1];
    std::memcpy(originalMatrix, currentMatrix, sizeof(int)*(N+1)*(N+1));
    this->originalVertices = new std::vector<seqGraph::Vertex*>();
    for (auto item : *this->graph->getVertices()) {
        this->originalVertices->push_back(new seqGraph::Vertex(*item));
    }
//    this->originalJunctions = new std::vector<seqGraph::Junction*>();
//    for (auto item : *this->graph->getJunctions()) {
//        this->originalJunctions->push_back(new seqGraph::Junction(*item));
//    }
//    this->originalVertices = this->graph->getVertices();
}

matching::~matching() {
//    delete this->matched;
    for (auto item: *originalVertices) {
        delete item;
    }
    originalVertices->clear();
    originalVertices->shrink_to_fit();
    for (int i = 0; i < N+1; ++i) {
        free(originalMatrix[i]);
    }
    free(originalMatrix);
}



void matching::resetGraph(seqGraph::Graph* g) {
    free(this->graph);
    std::cout<<"free done"<<std::endl;
    this->graph =  g;
    N = 2 *  g->getVCount();
    this->matched = new int[N + 1];
    std::fill_n(this->matched,N+1, -1);
//    for (int i = 0 ; i < N + 1; i++) this->matched[i] = -1;
    currentMatrix = this->graph->getConjugateMatrix();
    this->cyclePaths.clear();
    this->cyclePaths.shrink_to_fit();
}

bool matching::kmDfs(int u, bool visity[],bool visitx[], std::set<int>* pre, double ex[], double ey[], double slack[]) {
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
                pre->insert(u);
                pre->insert(conjugateIdx(i));
                matched[i] = u;
                matched[conjugateIdx(u)] = conjugateIdx(i);
                return true;
            }
        } else if(slack[i] > ex[u] + ey[i] - matrix[u][i]) {
//            auto tt1 = ex[u];
//            auto tt2 = ey[i];
//            auto tt3 = matrix[u][i];
//            if (ex[u] + ey[i]  - matrix[u][i] >= ZERO)
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

    std::fill_n(ex,N+1,0);
    std::fill_n(ey,N+1,0);
//    std::fill_n(ex,N+1,0);

//    memset(ey,0,sizeof ey);
//    memset(ex,0,sizeof ex);
//    memset(visity,0,sizeof visity);

    for( int i = 1 ; i < N +1 ; i++ ){
        ex[i] = *std::max_element(matrix[i], matrix[i]+N+1);
//        for( int j = 1 ; j <= N ; j++ ){
//            if( ex[i] < matrix[i][j] ) ex[i] = matrix[i][j];
//        }
    }
    auto* pre = new std::set<int>();
    for (int i = 1; i < N + 1; i++) {
        if (i == 18013) {
            int mmk = 1;
        }
        std::fill_n(slack, N+1, 10000000);
//        for( int l = 1 ; l <=N ; l++ ) slack[l] = 1000;
        if (pre->find(i) != pre->end()) {
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
//            for(int k = 0; k < N + 1; k++) {
//                std::cout<<matched[k]<<"\t";
//            }
            std::cout<<std::endl;
        }
    }
    free(pre);
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

std::string matching::idx2Str(int idx) {
    int now = idx;
    int vIdx = (now + 1) / 2;
    char dir = now % 2 == 0 ? '-':'+';
    auto idStr = (*this->originalVertices)[vIdx - 1]->getId();
    return idStr+dir;
}
std::map<int, std::vector<int>*>* matching::resolvePath(std::map<int, std::vector<int>*>* prevPaths) {
    auto matrix = graph->getConjugateMatrix();
    bool visited[N+1];
    memset(visited, 0, sizeof visited);
    auto* resolvedPath = new std::map<int, std::vector<int> *>();
    auto* resPath = new std::map<int, std::vector<int> *>();
//    auto checkedC = checkConjugateMatch();
//    if(checkedC!= 0)
//        std::cout<<"Conjugate not checked\n";
//    else
//        std::cout<<"check conjugate done\n";
    std::vector<int> cyclePaths;
    for(int i = 1; i < N + 1; i++) {
        if (visited[i]) continue;
        auto* currentPath = new std::vector<int>();
        if (matrix[matched[i]][i] == 0){
            if (matrix[matched[conjugateIdx(i)]][conjugateIdx(i)] == 0) {
                currentPath->push_back(i); resolvedPath->emplace(i,currentPath);
                visited[i] = true;
                visited[conjugateIdx(i)] = true;
            } else {
                visited[i] = true;
            }
            continue;
//            currentPath->push_back(i); resolvedPath->emplace(i,currentPath); continue;
        }
        int now = i;
        int vIdx = (now + 1) / 2;
        char dir = now % 2 == 0 ? '-':'+';
//        currentPath->push_back((*this->graph->getVertices())[vIdx - 1]->getId() + dir);
        currentPath->push_back(now);
//        if ((*this->graph->getVertices())[vIdx - 1]->getId() == "52_0")
//            int tm = 99;

//        std::cout<< (*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
        visited[now] = true;
        visited[conjugateIdx(now)] = true;
        bool currentInsert = true;
        bool isCycle = false;
        while (matrix[matched[now]][now] != double(0)) {
//            如果连上之前断裂的路径
            if (resolvedPath->find(matched[now]) != resolvedPath->end()) {
                auto oldPath = (*resolvedPath)[matched[now]];
//                +"_"+ std::to_string(now)
                for (auto it = currentPath->rbegin(); it != currentPath->rend(); it++){
                    oldPath->insert(oldPath->begin(), *it);
                }
                resolvedPath->erase(matched[now]);
//                auto extendPath = this->addPrevPath(prevPaths, oldPath);
                resolvedPath->emplace(i, oldPath);
                currentInsert = false;
                break;
            }
//            环状路径
            if(visited[matched[now]]){
                isCycle = true;
//                if (matched[now] == i) currentPath->push_back(-1);
                currentPath = breakCycle(currentPath);
                break;
            }
            visited[matched[now]] = true;
            visited[conjugateIdx(matched[now])] = true;
//            vIdx = (matched[now] + 1) / 2;
//            dir = matched[now] % 2 == 0 ? '-':'+';
//            std::cout<<(*this->graph->getVertices())[vIdx - 1]->getId()<<dir<<'\t';
            currentPath->push_back(matched[now]);
            now = matched[now];
        }
        if(currentInsert) {
            resolvedPath->emplace((*currentPath)[0], currentPath);
            if (isCycle) {
                this->cyclePaths.push_back((*currentPath)[0]);
            }
        }
    }
    for (auto path : *resolvedPath) {
        auto extendPath = this->addPrevPath(prevPaths, path.second);
        resPath->emplace(extendPath->front(),extendPath);
//        path.second->clear();
//        path.second->shrink_to_fit();
    }
    free(resolvedPath);
//    if(prevPaths != nullptr) {
//        for(auto item: *prevPaths) {
//            item.second->clear();
//            item.second->shrink_to_fit();
//        }
//    }
    free(prevPaths);

    free(matched);
    return resPath;
}

int conjugateIdx(int idx) {
    int cI;
    if(idx == 0) return idx;
    cI = idx%2 == 0? idx-1:idx+1;
    return cI;
}

double* mergePath(std::vector<int>* p1, std::vector<int>* p2, double** matrix, double* result) {
//    auto result = new double[4];
//    for(int i = 0 ; i< 4 ; i ++) result[i] = 0;

    for (auto ip1: *p1) {
        if (ip1 == -1) continue;
        for (auto ip2: *p2) {
            if (ip2 == -1) continue;
            result[0] += matrix[ip2][ip1];
            result[1] += matrix[conjugateIdx(ip2)][ip1];
            result[2] += matrix[ip2][conjugateIdx(ip1)];
            result[3] += matrix[conjugateIdx(ip2)][conjugateIdx(ip1)];
//            if (ip1 % 2 == 1) {
//                if (ip2 % 2 == 1) {//            ++
//                    result[0] += matrix[ip1][ip2];
//                } else {
//                    result[1] += matrix[ip1][ip2];
//                }
//            } else {
//                if (ip2 % 2 == 1) {//            -+
//                    result[2] += matrix[ip1][ip2];
//                } else {
//                    result[3] += matrix[ip1][ip2];
//                }
//            }
        }
    }
    return result;
}

void matching::reconstructMatrix(std::map<int, std::vector<int>*>* paths) {
    auto* resultG = new seqGraph::Graph();
//    auto tm = resultG->getConjugateMatrix() == nullptr;
    auto* values = new double[4];
    for (auto iPath: *paths) {
        resultG->addVertex(std::to_string(iPath.second->front()),"xx",1,2,1,1,2);
    }
    seqGraph::Vertex* v1;
    seqGraph::Vertex* v2;
//    for(auto junc : *originalJunctions) {
//        int i = junc->getSource()->getIdx();
//        int j = junc->getTarget()->getIdx();
//        int sDir = junc->getSourceDir();
//        int tDir = junc->getTargetDir();
//        double weightValue = junc->getWeight()->getCopyNum();
//        if (sDir == '+') {
//            if (tDir == '+') {
//                this->ConjugateMatrix[2*j + 1][2*i + 1] = weightValue;
//                this->ConjugateMatrix[2*(i+1)][2*(j+1)] = weightValue;
//            } else {
//                this->ConjugateMatrix[2*(j + 1)][2*i+1] = weightValue;
//                this->ConjugateMatrix[2*(i+1)][2*j+1] = weightValue;
//            }
//        } else {
//            if (tDir == '+') {
//                this->ConjugateMatrix[2*j+1][2*(i + 1)] = weightValue;
//                this->ConjugateMatrix[2*i+1][2*(j+1)] = weightValue;
//            } else {
//                this->ConjugateMatrix[2*i+1][2*j+1] = weightValue;
//                this->ConjugateMatrix[2*(j+1)][2*(i+1)] = weightValue;
//            }
//        }
//    }
    for (auto iPath: *paths) {
        v1 = resultG->getVertexByIdQ(std::to_string(iPath.second->front()));
        for (auto jPath: *paths) {
            if (iPath.first == jPath.first) continue;
            v2 = resultG->getVertexByIdQ(std::to_string(jPath.second->front()));
            if (iPath.second->size() == 1 && jPath.second->size() == 1) {
                int i = iPath.second->front();
                int j = jPath.second->front();
                auto matrix = this->originalMatrix;
                values[0] += matrix[j][i];
                values[1] += matrix[conjugateIdx(j)][i];
                values[2] += matrix[j][conjugateIdx(i)];
                values[3] += matrix[conjugateIdx(j)][conjugateIdx(i)];
            } else {
                mergePath(iPath.second, jPath.second, this->originalMatrix, values);
            }
            if (values[0] == 0 && values[1]==0 && values[2]==0 && values[3]==0) continue;
//            std::string v1Str, v2Str;
//            for (auto item : *iPath.second) {
//                v1Str+=idx2Str(item);
//            }
//            for (auto item : *jPath.second) {
//                v2Str+= idx2Str(item);
//            }

            resultG->addJunction(v1, v2, '+', '+', values[0], 1 , 1);
            resultG->addJunction(v1, v2, '+', '-', values[1], 1 , 1);
            resultG->addJunction(v1, v2, '-', '+', values[2], 1 , 1);
            resultG->addJunction(v1, v2, '-', '-', values[3], 1 , 1);
        }
    }
    delete[] values;
    std::cout<<"start reset"<<std::endl;
    resetGraph(resultG);
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

std::vector<int>* matching::breakCycle(std::vector<int> * cyclePath) {
    auto res = new std::vector<int>();
    auto matrix = this->getMatrix();
    int minIdx = cyclePath->size() - 1;
    double minW = INF;
    for (int i = 0; i < cyclePath->size(); i++) {
        if (i == cyclePath->size() - 1 && matrix[i][0] < minW) {
            minIdx = cyclePath->size() - 1;
        }
        if (matrix[i][i+1] < minW) {
            minW = matrix[i][i+1];
            minIdx = i;
        }
    }
    if (minIdx == cyclePath->size() - 1) return cyclePath;
    for (int i = minIdx + 1; i< cyclePath->size(); i++) {
        res->push_back((*cyclePath)[i]);
    }
    for (int i = 0; i <= minIdx; ++i) {
        res->push_back((*cyclePath)[i]);
    }
    return res;
}

std::vector<int>* matching::addPrevPath(std::map<int, std::vector<int>*>* prevPaths, std::vector<int>* curPath) {
    if (prevPaths == nullptr) return curPath;
    auto res = new std::vector<int>();
    for(auto item : *curPath) {
//        if (item == 0)
//            int tmi = 9;
//        if (item == -1) continue;
        int vIdx = (item + 1) / 2;
        int dir = item % 2;
        auto prevIdx = std::stoi((*this->graph->getVertices())[vIdx - 1]->getId());
//        if (prevIdx == 2980)
//            int mkh = 9;
        if(prevPaths->find(prevIdx) != prevPaths->end()) {
            auto pPath = (*prevPaths)[prevIdx];
//            if +
            if (dir == 1) {
                for (auto pItem: *pPath) {
//                    if (pItem == 0)
//                        int kh = 0;
                    res->push_back(pItem);
                }
            } else { // if -, reverse add
                for (auto it = pPath->rbegin(); it != pPath->rend(); ++it) {
//                    if (*it == -1)
//                        int kh = 0;
                    res->push_back(conjugateIdx(*it));
                }
            }
        } else {
            std::cout<<"error, prev path not found "<<prevIdx<<std::endl;
        }
    }
    return res;
}