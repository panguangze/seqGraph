//
// Created by caronkey on 10/5/2021.
//

#include "../include/Vertex.h"

using namespace seqGraph;

Vertex::Vertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double aCredibility, int aCopyNum) {
    this->Id = mId;
    this->orphan = true;
    this->chrom = aChrom;
    this->start = aStart;
    this->end = aEnd;
    this->credibility = aCredibility;
//    create weight
    weight = new Weight(aCoverage);
    weight->setCopyNum(aCopyNum);
//    create endpoints
    this->EP5 = new EndPoint(weight, _LEFT_TOP_, mId);
    this->rEP5 = new EndPoint(weight, _RIGHT_BOTTOM_, mId);
    this->EP3 = new EndPoint(weight, _RIGHT_TOP_, mId);
    this->rEP3 = new EndPoint(weight, _LEFT_BOTTOM_, mId);

    this->EP3->setVertex(this);
    this->rEP3->setVertex(this);
    this->EP5->setVertex(this);
    this->rEP5->setVertex(this);
    this->setMateEP();

    //    create inner edges, inner edges has same weight with vertex weigh
//    auto *topEdge = new Edge(this->EP5, this->EP3, aCoverage, _INNER_EDGE_);
//    auto *bottomEdge = new Edge(this->rEP5, this->rEP3, aCoverage, _INNER_EDGE_);
}


Vertex::~Vertex() {

}

bool Vertex::operator==(const Vertex &rhs) const {
    return this->Id == rhs.Id;
};

const std::string Vertex::getId() const {
    return Id;
}

void Vertex::setId(const int mId) {
    Vertex::Id = mId;
}

void Vertex::setIdx(const int idx) {
    this->idx = idx;
}

int Vertex::getIdx() const {
    return this->idx;
}
Weight *Vertex::getWeight() const {
    return weight;
}

void Vertex::setWeight(Weight *weight) {
    Vertex::weight = weight;
}

EndPoint *Vertex::getEp3() const {
    return EP3;
}

void Vertex::setEp3(EndPoint *ep3) {
    EP3 = ep3;
}

EndPoint *Vertex::getEp5() const {
    return EP5;
}

void Vertex::setEp5(EndPoint *ep5) {
    EP5 = ep5;
}

EndPoint *Vertex::getRep3() const {
    return rEP3;
}

void Vertex::setRep3(EndPoint *rEp3) {
    rEP3 = rEp3;
}

EndPoint *Vertex::getRep5() const {
    return rEP5;
}

void Vertex::setRep5(EndPoint *rEp5) {
    rEP5 = rEp5;
}


bool Vertex::isOrphan() const {
    return orphan;
}

void Vertex::setOrphan(bool orphan) {
    Vertex::orphan = orphan;
}

void Vertex::restoreCopy() {
    this->weight->restore();
}

void Vertex::backupCopy() {
    this->weight->backup();
}

void Vertex::setHasLowerBoundLimit() { hasLowerBoundLimit = true; }
void Vertex::resetHasLowerBoundLimit() { hasLowerBoundLimit = false; }
void Vertex::checkLowerBound() { hasLowerBoundLimit = this->weight->getCopyNum() >= 1; }

const std::string Vertex::getChrom() const {
    return chrom;
}

void Vertex::setChrom(const std::string mChrom) {
    Vertex::chrom = mChrom;
}

int Vertex::getStart() const {
    return start;
}

void Vertex::setStart(int mStart) {
    Vertex::start = mStart;
}

int Vertex::getEnd() const {
    return end;
}

void Vertex::setEnd(int mEnd) {
    Vertex::end = mEnd;
}

double Vertex::getCredibility() const {
    return credibility;
}

void Vertex::setCredibility(double mCredibility) {
    Vertex::credibility = mCredibility;
}

bool Vertex::isHasLowerBoundLimit() {
    return this->hasLowerBoundLimit;
}

bool Vertex::hasCopy() {
    return this->EP3->hasCopy();
}

double Vertex::getInCoverage(){
    return this->EP5->getInCoverage();
}
double Vertex::getOutCoverage(){
    return this->EP3->getOutCoverage();
}

void Vertex::setMateEP() {
    this->EP5->setMateEp(this->EP3);
    this->EP3->setMateEp(this->EP5);
    this->rEP5->setMateEp(this->rEP3);
    this->rEP3->setMateEp(this->rEP5);
}
