//
// Created by caronkey on 12/5/2021.
//

#ifndef SEQMAP_JUNCTION_H
#define SEQMAP_JUNCTION_H
#define _POSITIVE_DIR_ '+'
#define _NEGATIVE_DIR_ '+'

#include "Vertex.h"
#include "Edge.h"

namespace seqGraph {
    class Edge;
    class Vertex;
    class Junction {
    protected:
        char sourceDir;
        char targetDir;
        Vertex *source;
        Vertex *target;
        Weight *weight;
        bool hasLowerBoundLimit;

        Edge *oEdge; // original edge
        Edge *cEdge; // conjugate edge
        bool inferred;

    public:
        Junction(Vertex *sourceVertex, Vertex *targetVertex, char sourceDir, char targetDir, double copyNum,
                 double coverage, bool aIsBounded);

        ~Junction();

        void junctionToEdge();

        bool operator==(const Junction &otherJunc) const;

        char getSourceDir() const;

        char getTargetDir() const;

        Weight *getWeight() const;

        Edge *getOEdge() const;

        Edge *getCEdge() const;

        Vertex *getSource() const;

        Vertex *getTarget() const;

        bool isInferred() const;

        void setInferred(bool inferred);

        bool isHasLowerBoundLimit() const;

        void setHasLowerBoundLimit();
        void resetHasLowerBoundLimit();
        void checkLowerBound();

        std::string getInfo();


        void restoreCopy();
        void backupCopy();

        bool hasCopy();
    };
}


#endif //SEQMAP_JUNCTION_H
