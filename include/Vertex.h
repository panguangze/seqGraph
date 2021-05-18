//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_VERTEX_H
#define SEQMAP_VERTEX_H


#include <string>
#include "EndPoint.h"
#include "Weight.h"

namespace seqGraph {
    class EndPoint;
    class Vertex {
    protected:
        std::string Id;
        std::string chrom;
        int start;
        int end;
        double credibility;

        Weight *weight;
        EndPoint *EP3;
        EndPoint *EP5;
        EndPoint *rEP3;
        EndPoint *rEP5;
        bool orphan;
        bool hasLowerBoundLimit;
    public:
        Vertex(std::string mId, std::string aChrom, int aStart, int aEnd,double aCoverage, double mCredibility, int aCopyNum);

        ~Vertex();

        const std::string getId() const;

        void setId(const std::string mId);

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

        bool operator==(const Vertex &rhs) const;

        void restoreCopy();
        void backupCopy();

        void setHasLowerBoundLimit();
        bool isHasLowerBoundLimit();
        void resetHasLowerBoundLimit();
        void checkLowerBound();
        bool hasCopy();

        double getInCoverage();
        double getOutCoverage();
    };
}


#endif //SEQMAP_VERTEX_H
