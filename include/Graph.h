//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_GRAPH_H
#define SEQMAP_GRAPH_H

#include "Vertex.h"
#include "Junction.h"
#include  <set>
#include "map"
extern std::string SUB_ONLY;

namespace seqGraph {
    typedef std::vector<EndPoint *> EndPointPath;


    //    CRS for sparse matrix
    class SparseMatrix {
    public:
        int N;

        std::vector<float> values;
        std::vector<int> IA;
        std::vector<int> JA;

        explicit SparseMatrix ();
        ~SparseMatrix();

        float getIJ(int i, int j);
//        void addValue(int i, int j, float v);
        void debugPrint();
        float getIRowMax(int i);
    };


    class Graph {
    protected:
        float mAvgCoverage;
        std::vector<Vertex *> *vertices;
        std::vector<Junction *> *junctions;
        Vertex *source;
        Vertex *sink;
        float ** ConjugateMatrix;
        std::map<std::string, int>* verticesIdx;
        std::map<std::string, int>* junctionIdx;
        std::vector<std::set<int>* > connectedJunctionsIdx;
        SparseMatrix sparseMatrix;

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

        Vertex *addVertex(std::string mId, std::string aChrom, int aStart, int aEnd,float aCoverage, float mCredibility, int aCopyNum);

        Junction *
        addJunction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, float copyNum,
                    float coverage, bool aIsBounded);

        Junction *
        addJunction(std::string& sourceId, std::string& targetId, char sourceDir, char targetDir, float copyNum,
                    float coverage, bool aIsBounded);

        Junction * addJunction(EndPoint* ep3, EndPoint* ep5, float copyNum, float converage, bool isBounded);

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

        float getMAvgCoverage() const;

        void setMAvgCoverage(float mAvgCoverage);

        Vertex *getSource() const;

        void setSource(std::string& Id);

        Vertex *getSink() const;

        void setSink(std::string& Id);

        std::vector<Vertex *> *getVertices() const;

        std::vector<Junction *> *getJunctions() const;

//        float ** getConjugateMatrix();

        SparseMatrix& getConjugateMatrix();
    };
}

#endif //SEQMAP_GRAPH_H
