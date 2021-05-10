//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_VERTEX_H
#define SEQMAP_VERTEX_H


#include <string>
#include "VertexMeta.h"
#include "EndPoint3.h"
#include "EndPoint5.h"
#include "Weight.h"

class Vertex {
protected:
    int mId;    // segment id
    std::string mChrom;
    int mStart;
    int mEnd;

    double mCredibility;
    bool mHasLowerBoundLimit;  // used in ILP processing, if true then lower bound is 1, otherwise 0
    bool mIsOrphan;
    bool mHasCheckedOrphan;

    Weight * mWeight;

    EndPoint3* EP3;
    EndPoint5* EP5;
    EndPoint3* rEP3;
    EndPoint5* rEP5;

};


#endif //SEQMAP_VERTEX_H
