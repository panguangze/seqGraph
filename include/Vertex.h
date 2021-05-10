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
    std::string mId;    // segment id

    EndPoint3* EP3;
    EndPoint5* EP5;
    EndPoint3* rEP3;
    EndPoint5* rEP5;

};


#endif //SEQMAP_VERTEX_H
