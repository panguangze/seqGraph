//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_VERTEX_H
#define SEQMAP_VERTEX_H


#include <string>
#include "EndPoint.h"
#include "Weight.h"
#include <vector>
#include "Junction.h"

namespace seqGraph {
    class EndPoint;
    class Junction;
    class Vertex {
    protected:
        std::string Id;
        std::string chrom;
        int start;
        int end;
        int length;
        double depth;
        double credibility;
        int idx;
        int copy_idx; //

        bool visited;
    public:
        bool isVisited() const;

        void setIsVisited(bool isVisited);

    protected:

        Weight *weight;
        EndPoint *EP3;
        EndPoint *EP5;
        EndPoint *rEP3;
        EndPoint *rEP5;
        bool orphan;
        bool hasLowerBoundLimit;
        std::vector<Junction*> nextJuncs;
        std::vector<Junction*> prevJuncs;
    public:
        const std::vector<Junction *> &getNextJuncs() const;

        const std::vector<Junction *> &getPrevJuncs() const;
    public:
        Vertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum);
        bool sameVertex(const seqGraph::Vertex & v2);

        ~Vertex();

        void setNextJunc(Junction*);
        void setPrevJunc(Junction*);


        const std::string getId() const;

        bool operator ==(const Vertex &otherVertex) const;

        bool operator >(const Vertex &otherVertex) const;
        bool operator <(const Vertex &rhs) const;


        void setId(const int mId);

        void setIdx(const int idx);

        int getIdx() const;

        Weight *getWeight() const;

        void setWeight(Weight *weight);

        EndPoint *getEp3() const;

        void setEp3(EndPoint *ep3);

        EndPoint *getEp5() const;

        void setEp5(EndPoint *ep5);

        EndPoint *getRep3() const;

        void setRep3(EndPoint *rEp3);

        EndPoint *getRep5() const;

        void setRep5(EndPoint *rEp5);

        bool isOrphan() const;

        void setOrphan(bool orphan);

        const std::string getChrom() const;

        void setChrom(const std::string mChrom);

        int getStart() const;

        void setStart(int mStart);

        int getEnd() const;

        void setEnd(int mEnd);

        double getCredibility() const;

        void setCredibility(double mCredibility);


        void restoreCopy();
        void backupCopy();

        void setHasLowerBoundLimit();
        bool isHasLowerBoundLimit();
        void resetHasLowerBoundLimit();
        void checkLowerBound();
        bool hasCopy();

        double getInCoverage();
        double getOutCoverage();

        void setMateEP();
    };
}


#endif //SEQMAP_VERTEX_H
