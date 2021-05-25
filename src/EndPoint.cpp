//
// Created by caronkey on 10/5/2021.
//

#include "../include/EndPoint.h"

using namespace seqGraph;

EndPoint::EndPoint(Weight *weight, int type, int vId) {
    this->weight = weight;
    this->type = type;
    this->visited = false;
    this->shortestPrevEdge = nullptr;
    this->edges = new std::vector<Edge *>();
    this->vId = vId;
}

EndPoint::~EndPoint() {

}

void EndPoint::addEdge(Edge *edge) {
    this->edges->push_back(edge);
}

int EndPoint::getType() const {
    return type;
}

void EndPoint::setType(int type) {
    EndPoint::type = type;
}

std::string EndPoint::getInfo() {
    return std::to_string(this->vId) + std::to_string(this->type);
}

Edge *EndPoint::getOneNextEdge(bool isTraversing) {
//    if EP 5 get its mateEP
    EndPoint *currentEP;
    if (this->type == _LEFT_TOP_ || this->type == _RIGHT_BOTTOM_) currentEP = this->mateEP;
    else currentEP = this;
    for (Edge *e : *(currentEP->edges)) {
        if (!e->isVisited()) {
            if (isTraversing) {
                if (e->hasCopy()) {
                    return e;
                } else {
                    continue;
                }
            } else {
                return e;
            }
        }
    }
    return nullptr;
}

bool EndPoint::hasCopy() {
    return this->weight->getCopyNum() >= 1;
}

void EndPoint::traverse() {
    this->weight->decreaseCopyNum();
}

bool EndPoint::isVisited() const {
    return visited;
}

void EndPoint::setVisited(bool visited) {
//    this and its mate visited
    EndPoint::visited = visited;
    this->mateEP->visited = visited;
}

Weight *EndPoint::getWeight() const {
    return weight;
}

Edge *EndPoint::getShortestPrevEdge() const {
    return shortestPrevEdge;
}

void EndPoint::setShortestPrevEdge(Edge *shortestPrevEdge) {
    EndPoint::shortestPrevEdge = shortestPrevEdge;
    this->mateEP->shortestPrevEdge = shortestPrevEdge;
}

const int EndPoint::getVId() const {
    return vId;
}

EndPoint *EndPoint::getMateEp() const {
    return mateEP;
}

std::vector<Edge *> *EndPoint::getInEdges() {
    if (this->type == _LEFT_TOP_ || this->type == _RIGHT_BOTTOM_) return this->edges;
    if (this->type == _LEFT_BOTTOM_ || this->type == _RIGHT_TOP_) return this->mateEP->edges;
}

std::vector<Edge *> *EndPoint::getOutEdges() {
    if (this->type == _LEFT_TOP_ || this->type == _RIGHT_BOTTOM_) return this->mateEP->edges;
    if (this->type == _LEFT_BOTTOM_ || this->type == _RIGHT_TOP_) return this->edges;
}

double EndPoint::getInCoverage() {
    double inCoverage = 0;
    for (Edge * e : *(this->getInEdges())) {
        inCoverage += e->getWeight()->getCoverage();
    }
    return inCoverage;
}
double EndPoint::getOutCoverage() {
    double outCoverage = 0;
    for (Edge * e : *(this->getOutEdges())) {
        outCoverage += e->getWeight()->getCoverage();
    }
    return outCoverage;
}

Vertex *EndPoint::getVertex() const {
    return vertex;
}

void EndPoint::setVertex(Vertex *vertex) {
    EndPoint::vertex = vertex;
}

void EndPoint::setMateEp(EndPoint *mateEp) {
    mateEP = mateEp;
}
