//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_GRAPH_H
#define SEQMAP_GRAPH_H

#include "Vertex.h"
#include "Junction.h"
#include  <set>
#include "map"

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
        std::map<std::string, int>* verticesIdx;
        std::map<std::string, int>* junctionIdx;
        std::vector<std::set<int>* > connectedJunctionsIdx;
    public:
        bool isReconstructed;
        Graph();

        ~Graph();
        Graph* getSubgraph(int i);
        void originalGraphFree();
        inline int getVCount() const {
            return this->vertices->size();
        }

        inline int subGraphCount() {
            return this->connectedJunctionsIdx.size();
        }

        inline int getJuncSize() {
            return this->junctionIdx->size();
        }

        Vertex *addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum);

        Junction *
        addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction *
        addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, double copyNum,
                    double coverage, bool aIsBounded);

        Junction * addJunction(EndPoint* ep3, EndPoint* ep5, double copyNum, double converage, bool isBounded);

        Junction *
        addJunction(Junction* j);

        Vertex *addVertex(Vertex *v);

        Vertex *getVertexById(std::string Id);
        Vertex *getVertexByIdQ(std::string Id);
        Vertex *getVertexByIdx(int idx);

        bool doesPathExists(EndPoint *sourceEndPoint, EndPoint *sinkEndpoint);

        int BFS_EndPoint(EndPoint *sourceEndpoint, EndPoint *sinkEndpoint);

        void BFS_Vertices(Vertex* vertex,std::set<int>* connectedIdx);

        Junction* doesJunctionExist(Junction *junction);

        Junction* doesJunctionExist(Vertex& v1, Vertex& v2, char v1d, char v2d);

        void resetJunctionVisitFlag();

        void parseConnectedComponents();

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
