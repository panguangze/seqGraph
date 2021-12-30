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
        double ** ConjugateMatrix;
    public:
        Graph();

        ~Graph();

        inline int getVCount() const {
            return this->vertices->size();
        }

        Vertex *addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum);

        Junction *
        addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction *
        addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction * addJunction(EndPoint* ep3, EndPoint* ep5, double copyNum, double converage, bool isBounded);

        Vertex *getVertexById(std::string Id);

        bool doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint);

        int BFS(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint);

        bool doesJunctionExist(Junction *junction);

        void resetJunctionVisitFlag();

        void resetVertexVisitFlag();

        void resetShortestPrevEdge();

        double getMAvgCoverage() const;

        void setMAvgCoverage(double mAvgCoverage);

        Vertex *getSource() const;

        void setSource(std::string& Id);

        Vertex *getSink() const;

        void setSink(std::string& Id);

        std::vector<Vertex *> *getVertices() const;

        std::vector<Junction *> *getJunctions() const;

        double ** getConjugateMatrix();
    };
}

#endif //SEQMAP_GRAPH_H
