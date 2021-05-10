//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_ENDPOINT3_H
#define SEQMAP_ENDPOINT3_H

#include <vector>
#include "EndPoint5.h"


class EndPoint3 {
protected:
    int *copyNumber;
    std::vector<EndPoint5 *> nextPoints;
};


#endif //SEQMAP_ENDPOINT3_H
