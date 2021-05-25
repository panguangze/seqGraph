//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_GRAPH_H
#define SEQMAP_GRAPH_H

#include "Vertex.h"
#include "Junction.h"

namespace seqGraph {
    typedef std::vector<EndPoint *> EndPointPath;

    class Graph {
    protected:
        double mAvgCoverage;
        std::vector<Vertex *> *vertices;
        std::vector<Junction *> *junctions;
        Vertex *source;
        Vertex *sink;
    public:
        Graph();

        ~Graph();

        Vertex *addVertex(int mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum);

        Junction *
        addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction *
        addJunction(int sourceId, int targetId, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction * addJunction(EndPoint* ep3, EndPoint* ep5, double copyNum, double converage, bool isBounded);

        Vertex *getVertexById(int Id);

        bool doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint);

        int BFS(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint);

        bool doesJunctionExist(Junction *junction);

        void resetJunctionVisitFlag();

        void resetVertexVisitFlag();

        void resetShortestPrevEdge();

        double getMAvgCoverage() const;

        void setMAvgCoverage(double mAvgCoverage);

        Vertex *getSource() const;

        void setSource(int Id);

        Vertex *getSink() const;

        void setSink(int Id);

        std::vector<Vertex *> *getVertices() const;

        std::vector<Junction *> *getJunctions() const;
    };
}

#endif //SEQMAP_GRAPH_H
