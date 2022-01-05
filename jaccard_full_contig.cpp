//
// Created by caronkey on 27/11/2020.
//

#include "jaccard_full_contig.h"

// arg1 输入bam，arg2输出jaccard
int main(int argc, char **argv) {
    htsFile *in;

    if ((in = hts_open(argv[1], "r")) == nullptr) {
        fprintf(stderr, "Error opening '%s'\n", argv[1]);
        return -1;
    }
    readBAM(in, argv[2], 100);
}

void initIMap (sam_hdr_t *hdr,Interactions& iMap) {
    int nRef = sam_hdr_nref(hdr);
    std::string refName, revRef;
    for(int i = 0; i< nRef; i++) {
        refName = sam_hdr_tid2name(hdr, i);
        revRef = refName+" -";
        std::unordered_map<std::string, int> tmp;
        std::unordered_map<std::string, int> rvtmp;
        iMap[refName] = tmp;
        iMap[revRef]= rvtmp;
    }
}

void readBAM(htsFile *in, const char* out_file, int readsLen) {
    Interactions iMap;
    sam_hdr_t *hdr;
    bam1_t *b;
    int ret;
    if ((hdr = sam_hdr_read(in)) == nullptr) {
        fprintf(stderr, "[E::%s] couldn't read file header \n", __func__);
        return;
    }
    if ((b = bam_init1()) == nullptr) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    std::string readName, refName, mRefName;
    initIMap(hdr, iMap);
    int pos, mpos, refLen, mRefLen;
    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        if (ret < -1) {
            fprintf(stderr, "[E::%s] Error parsing input.\n", __func__);
            if (b) bam_destroy1(b);
            if (hdr) sam_hdr_destroy(hdr);
        }
        if (b->core.tid == -1 || b->core.mtid == -1)
            continue;
        refName = std::string(sam_hdr_tid2name(hdr, b->core.tid));
        mRefName = std::string(sam_hdr_tid2name(hdr, b->core.mtid));
        pos = b->core.pos;
        mpos = b->core.mpos;
        auto flags = b->core.flag;
        auto rev = flags & 0x10;
        auto mrev = flags & 0x20;
//        refLen = sam_hdr_tid2len(hdr, b->core.tid);
//        mRefLen = sam_hdr_tid2len(hdr, b->core.mtid);
        if (refName == mRefName)continue;
        if (rev)
            refName = refName.append(" -");
        if (!mrev)
            mRefName = mRefName.append(" -");
//        if (refLen <= readsLen) {
//            if (pos == 0)
//                refName = refName.append(" -");
//        } else if (pos < refLen/2 - readsLen/2) {
//            refName = refName.append(" -");
//        }
//        if (mRefLen <= readsLen) {
//            if (pos != 0)
//                mRefName = mRefName.append(" -");
//        } else if (pos >= mRefLen/2 - readsLen/2) {
//            mRefName = mRefName.append(" -");
//        }
        auto& refMap = iMap[refName];
        if (refMap.find(mRefName) == refMap.end()) {
            refMap[mRefName] = 1;
        } else {
            refMap[mRefName]+=1;
        }
    }
    std::ofstream fout(out_file);
    std::string prev;
    std::string next;
    for (auto& it: iMap) {
        for(auto& it2 : it.second) {
            if (it.first[-1] != '-') prev = it.first + " +";
            else prev = it.first;
            if (it2.first[-1] != '-') next = it2.first + " +";
            else next = it2.first;
            fout<<prev<<" "<<next<<" "<<it2.second<<"\n";
        }
//        if (it.second > 50) fout<<it.first.first<<"\t"<<it.first.second<<"\t"<<it.second<<"\t"<<"\n";
    }
    fout.close();
}

