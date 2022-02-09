//
// Created by caronkey on 10/5/2021.
//

#include "../include/Graph.h"
#include "../include/Exceptions.h"
#include <iostream>
#include <algorithm>
#include <queue>

using namespace seqGraph;

Graph::Graph() {
    this->mAvgCoverage = 0;
    this->vertices = new std::vector<Vertex *>();
    this->junctions = new std::vector<Junction *>();
    verticesIdx = new std::map<std::string, int>();
    junctionIdx = new std::map<std::string, int>();
    ConjugateMatrix = nullptr;
    source = nullptr;
    sink = nullptr;
}

Graph::~Graph() {
    for (auto item : *vertices) {
        delete item;
    }
    for (auto item : *junctions) {
        delete item;
    }
    vertices->clear();
    vertices->shrink_to_fit();
    junctions->clear();
    junctions->shrink_to_fit();
    verticesIdx->clear();
    junctionIdx->clear();
//    if (ConjugateMatrix != nullptr)
}

Graph* Graph::getSubgraph(int i) {
    auto subG = new Graph();
    auto idxs = connectedJunctionsIdx[i];
    for(auto j : *idxs) {
        auto junc = (*this->junctions)[j];
        subG->addJunction(junc);
    }
    return subG;
}
Vertex *Graph::addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum) {
//    create vertex add push
    auto v1 = this->getVertexById(mId);
    if(v1 != nullptr) return v1;
    auto *vertex = new Vertex(mId, aChrom, aStart, aEnd, aCoverage, mCredibility, aCopyNum);
    vertex->setIdx(this->vertices->size());
    this->vertices->push_back(vertex);
    this->verticesIdx->emplace(mId, vertex->getIdx());
    return vertex;
}

Junction *Graph::addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                             double coverage, bool aIsBounded) {
//    auto *junction = new Junction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
    auto jun = this->doesJunctionExist(*sourceVertex, *targetVertex, sourceDir, targetDir);
    if (jun == nullptr) {
        auto  k = sourceVertex->getId()+ targetVertex->getId()+sourceDir+targetDir;
//        throw DuplicateJunctionException(junction);
        auto *junction = new Junction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
        junction->junctionToEdge();
        junction->setIdx(this->junctions->size());

        this->junctions->push_back(junction);
        this->junctionIdx->emplace(k, junctions->size());
    //    set vertex not orphan
        sourceVertex->setOrphan(false);
        targetVertex->setOrphan(false);
        sourceVertex->setNextJunc(junction);
        targetVertex->setPrevJunc(junction);
        return junction;
    }
    return jun;
}

Junction *
Graph::addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, double copyNum,
                   double coverage ,bool aIsBounded) {
    Vertex *sourceVertex = this->getVertexById(sourceId);
    Vertex *targetVertex = this->getVertexById(targetId);
    return this->addJunction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
}

Junction * Graph::addJunction(EndPoint *ep3, EndPoint *ep5, double copyNum, double coverage, bool isBounded) {
    Vertex *sourceVertex = ep3->getVertex();
    Vertex *targetVertex = ep5->getVertex();

    char sourceDir = ep3->getType() == _RIGHT_TOP_ ? '+' : '-';
    char targetDir = ep3->getType() == _LEFT_TOP_ ? '+' : '-';
    return this->addJunction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, isBounded);
}

Vertex *Graph::getVertexById(std::string Id) {
    if (verticesIdx->find(Id) != verticesIdx->end()){
        auto idx = (*verticesIdx)[Id];
        return (*this->vertices)[idx];
    }
//    for (auto *vertex : *this->vertices) {
//        if (vertex->getId() == Id) return vertex;
//    }
//    throw VertexDoesNotExistException(Id);
    return nullptr;
}


int Graph::BFS_EndPoint(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint) {
    std::queue<EndPoint *> EPQueue;
    EPQueue.push(sourceEndpoint);

    while (!EPQueue.empty()) {
        EndPoint *currentEP = EPQueue.front();
        EPQueue.pop();

        for (Edge *e : *(currentEP->getOutEdges())) {
            // if (e->hasCopy()) {
            EndPoint *nextEP = e->getTarget();
            if (!nextEP->isVisited()) {
                nextEP->setVisited(true);
                nextEP->setShortestPrevEdge(e);
                EPQueue.push(nextEP);
            }
        }
    }
    this->resetVertexVisitFlag();

    Edge *prevEdge;
    EndPoint *prevEP;
    EndPoint *currentEP = sinkEndpoint;
    int found;
    while (true) {
        prevEdge = currentEP->getShortestPrevEdge();
        if (prevEdge == nullptr) {
            found = -1;
            break;
            // return -1;  // no path from aStartVertex to aTargetVertex
        }

        prevEP = prevEdge->getSource();
        if (prevEP == sourceEndpoint) {
            found = 0;
            break;
            // return 0;
        }
        currentEP = prevEP;
    }
    this->resetShortestPrevEdge();
    return found;
}

void Graph::BFS_Vertices(Vertex* vertex, std::vector<int>* connectedIdx){
    std::queue<Vertex*> VQueuen;
    VQueuen.push(vertex);

    while (!VQueuen.empty()) {
        auto currentV = VQueuen.front();
        VQueuen.pop();
        for (auto e : currentV->getNextJuncs()) {
            // if (e->hasCopy()) {
            auto *nextV = e->getTarget();
            if (!nextV->isVisited()) {
                connectedIdx->push_back(e->getIdx());
                nextV->setIsVisited(true);
                VQueuen.push(nextV);
            }
        }
        for (auto e : currentV->getPrevJuncs()) {
            // if (e->hasCopy()) {
            auto *nextV = e->getSource();
            if (!nextV->isVisited()) {
                connectedIdx->push_back(e->getIdx());
                nextV->setIsVisited(true);
                VQueuen.push(nextV);
            }
        }
    }
//    this->resetVertexVisitFlag();
}



Junction* Graph::doesJunctionExist(Junction *aJunction) {
    for(auto item: *this->junctions){
        if(*item == *aJunction) {
            return item;
        }
    }
    return nullptr;
}


Junction* Graph::doesJunctionExist(Vertex& v1, Vertex& v2, char v1d, char v2d) {
    auto  k = v1.getId()+ v2.getId()+v1d+v2d;
    if (junctionIdx->find(k) != junctionIdx->end()) {
        auto idx = (*junctionIdx)[k];
        return (*junctions)[idx];
    }
//    for(auto item: *this->junctions){
//        auto& s1 = *item->getSource();
//        auto& s2 = *item->getTarget();
//        auto s1d = item->getSourceDir();
//        auto s2d = item->getTargetDir();
//        if ((s1 == v1 && s2 == v2 && v1d == s1d && v2d == s2d) || (s1 == v2 && s2 == v1 && v1d != s1d && v2d != s2d) )
//            return item;
//    }
    return nullptr;
}
bool Graph::doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint) {
    bool isReach = false;
    EndPoint *startEP = sourceEndPoint;
    EndPointPath EPStack;
    EPStack.push_back(startEP);
    Edge *nextEdge = startEP->getOneNextEdge();
    while (true) {
        if (nextEdge == NULL) {
            EPStack.pop_back();
            if (EPStack.empty()) {
                break;
            } else {
                startEP = EPStack.back();
                nextEdge = startEP->getOneNextEdge();
            }
        } else {
            EndPoint *nextEP = nextEdge->getTarget();
            if (nextEP == sinkEndpoint || nextEP == sinkEndpoint->getMateEp()) {
                isReach = true;
                break;
            }
            nextEdge->setVisited(true);
            EPStack.push_back(nextEP);
            startEP = nextEP;
            nextEdge = startEP->getOneNextEdge();
        }
    }
    this->resetJunctionVisitFlag();
    return isReach;
}

void Graph::resetJunctionVisitFlag() {
    for (Junction *junction : *this->junctions) {
        junction->getCEdge()->setVisited(false);
        junction->getOEdge()->setVisited(false);
    }
}

void Graph::resetVertexVisitFlag() {
    for (auto *vertex : *this->vertices) {
        vertex->getEp3()->setVisited(false);
        vertex->getRep3()->setVisited(false);
    }
}

void Graph::resetShortestPrevEdge() {
    for (auto *vertex : *this->vertices) {
        vertex->getEp3()->setShortestPrevEdge(nullptr);
        vertex->getRep3()->setShortestPrevEdge(nullptr);
    }
}

double Graph::getMAvgCoverage() const {
    return mAvgCoverage;
}

void Graph::setMAvgCoverage(double mAvgCoverage) {
    Graph::mAvgCoverage = mAvgCoverage;
}

Vertex *Graph::getSource() const {
    return source;
}

void Graph::setSource(std::string& Id) {
    this->source = this->getVertexById(Id);
}

Vertex *Graph::getSink() const {
    return sink;
}

void Graph::setSink(std::string& Id) {
    this->sink = this->getVertexById(Id);
}

std::vector<Vertex *> *Graph::getVertices() const {
    return vertices;
}

std::vector<Junction *> *Graph::getJunctions() const {
    return junctions;
}

double ** Graph::getConjugateMatrix(){
    if(this->ConjugateMatrix != nullptr) return this->ConjugateMatrix;
    int n = this->getVCount();
    this->ConjugateMatrix = new double *[2*n + 1];
    for(int i = 0; i < 2*n + 1; i++) {

        auto t = new double[2*n + 1];
        std::fill_n(t,(2*n + 1), 0);
        this->ConjugateMatrix[i] = t;
    }
    double initV = 0;
    for(auto junc : *junctions) {
        int i = junc->getSource()->getIdx();
        int j = junc->getTarget()->getIdx();
        int sDir = junc->getSourceDir();
        int tDir = junc->getTargetDir();
        double weightValue = junc->getWeight()->getCopyNum();
        if (sDir == '+') {
            if (tDir == '+') {
                this->ConjugateMatrix[2*j + 1][2*i + 1] = weightValue;
                this->ConjugateMatrix[2*(i+1)][2*(j+1)] = weightValue;
            } else {
                this->ConjugateMatrix[2*(j + 1)][2*i+1] = weightValue;
                this->ConjugateMatrix[2*(i+1)][2*j+1] = weightValue;
            }
        } else {
            if (tDir == '+') {
                this->ConjugateMatrix[2*j+1][2*(i + 1)] = weightValue;
                this->ConjugateMatrix[2*i+1][2*(j+1)] = weightValue;
            } else {
                this->ConjugateMatrix[2*i+1][2*j+1] = weightValue;
                this->ConjugateMatrix[2*(j+1)][2*(i+1)] = weightValue;
            }
        }
    }
    return ConjugateMatrix;
}

void Graph::parseConnectedComponents() {
    int i = 0;
    for(auto & vertice : *vertices) {
        if (i == vertices->size()) break;
        if (!vertice->isVisited()) {
            auto * connectedIdx = new std::vector<int>;
            BFS_Vertices(vertice, connectedIdx);
            i += connectedIdx->size();
            this->connectedJunctionsIdx.push_back(connectedIdx);
        }
    }
}

Vertex *Graph::addVertex(Vertex *v) {
    auto v1 = this->getVertexById(v->getId());
    if(v1 != nullptr) return v1;
    v->setIdx(this->vertices->size());
    this->vertices->push_back(v);
    this->verticesIdx->emplace(v->getId(), v->getIdx());
    return v;
}

Junction *Graph::addJunction(Junction *j) {
    this->addVertex(j->getSource());
    this->addVertex(j->getTarget());
    auto k = j->getSource()->getId()+ j->getTarget()->getId()+j->getSourceDir()+j->getTargetDir();
    this->junctions->push_back(j);
    j->setIdx(this->junctions->size());
    this->junctionIdx->emplace(k, junctions->size());
    return j;
}