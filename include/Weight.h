//
// Created by caronkey on 10/5/2021.
//

#ifndef SEQMAP_WEIGHT_H
#define SEQMAP_WEIGHT_H


class Weight {
protected:
    double mCoverage;           // coverage as adjusted by LP
    double mCopyNum;             // the copy number, coverage / average coverage of one copy
    double mCopyNumOriginal;
    double mCopyNumBackup;

    bool mIsInferred;

public:
    Weight(double aCoverage);
    ~Weight();

    double getCoverage();
    double getCopyNum();
    double getCopyNumBackup();
    void setCoverage(double aCoverage);
    void setCopyNum(double aCopyNum);
    void backup();
    void restore();

    void increaseCopyNum(double aIncrement = 1);
    void decreaseCopyNum(double aDecrement = 1);

    bool isInferred();
    void setInferred();
    void resetInferred();

    void print();
};


#endif //SEQMAP_WEIGHT_H
