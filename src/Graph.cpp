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
    source = nullptr;
    sink = nullptr;
}

Graph::~Graph() {

}

Vertex *Graph::addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum) {
//    create vertex add push

    auto *vertex = new Vertex(mId, aChrom, aStart, aEnd, aCoverage, mCredibility, aCopyNum);
    vertex->setIdx(this->vertices->size());
    this->vertices->push_back(vertex);
    return vertex;
}

Junction *Graph::addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                             double coverage, bool aIsBounded) {
    auto *junction = new Junction(sourceVertex, targetVertex, sourceDir, targetDir, copyNum, coverage, aIsBounded);
    if (this->doesJunctionExist(junction)) {
        throw DuplicateJunctionException(junction);
    }
    junction->junctionToEdge();
    this->junctions->push_back(junction);
//    set vertex not orphan
    sourceVertex->setOrphan(false);
    targetVertex->setOrphan(false);
    return junction;
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
    for (auto *vertex : *this->vertices) {
        if (vertex->getId() == Id) return vertex;
    }
//    throw VertexDoesNotExistException(Id);
    return nullptr;
}


int Graph::BFS(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint) {
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


bool Graph::doesJunctionExist(Junction *aJunction) {
    return std::find(this->junctions->begin(), this->junctions->end(), aJunction) != this->junctions->end();
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
        this->ConjugateMatrix[i] = new double[2*n + 1];
    }
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