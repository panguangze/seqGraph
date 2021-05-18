//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_ENDPOINT_H
#define SEQMAP_ENDPOINT_H
#define _LEFT_TOP_ 0
#define _RIGHT_TOP_ 1
#define _RIGHT_BOTTOM_ 2
#define _LEFT_BOTTOM_ 3

#include "Weight.h"
#include "Edge.h"
#include "Vertex.h"
#include <vector>
#include <string>

namespace seqGraph {
    class Edge;
    class Vertex;
    class EndPoint {
    protected:
        Weight *weight;
        EndPoint *mateEP;
        Vertex* vertex;
        std::string vId;
        std::vector<Edge *> *edges;
        int type; // EP type,0: left top, 1: right top, 2: right bottom,3: left bottom
        bool visited;
        Edge *shortestPrevEdge;
    public:
        EndPoint(Weight *weight, int type, std::string vId);

        ~EndPoint();

        void addEdge(Edge *edge);

        std::string getInfo();

        Edge *getOneNextEdge(bool isTraversing = false);

        bool hasCopy();

        void traverse();

        int getType() const;

        void setType(int type);

        bool isVisited() const;

        Weight *getWeight() const;

        void setVisited(bool visited);

        Edge *getShortestPrevEdge() const;

        void setShortestPrevEdge(Edge *shortestPrevEdge);

        const std::string getVId() const;

        EndPoint *getMateEp() const;

        std::vector<Edge *> *getInEdges();

        std::vector<Edge *> *getOutEdges();

        Vertex *getVertex() const;

        void setVertex(Vertex *vertex);

        double getInCoverage();
        double getOutCoverage();

    };
}


#endif //SEQMAP_ENDPOINT_H
