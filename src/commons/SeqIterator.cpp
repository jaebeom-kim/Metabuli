#include "SeqIterator.h"

SeqIterator::~SeqIterator() {
}
SeqIterator::SeqIterator(const LocalParameters &par) {

}

string SeqIterator::reverseComplement(string &read) const {
    int len = read.length();
    string out;
    for (int i = 0; i < len; i++) {
        out.push_back(iRCT[read[i]]);
    }
    reverse(out.begin(), out.end());
    return out;
}

char *SeqIterator::reverseComplement(char *read, size_t length) const {
    char *revCom = (char *) malloc(sizeof(char) * (length + 1));
    for (size_t i = 0; i < length; i++) {
        revCom[length - i - 1] = iRCT[read[i]];
    }
    revCom[length] = '\0';
    return revCom;
}

bool SeqIterator::compareMinHashList(priority_queue <uint64_t> list1, priority_queue <uint64_t> &list2, size_t length1,
                                     size_t length2) {
    float lengthRatio = float(length2) / float(length1);
    float identicalCount = 0;
    float list1Size = list1.size();
    while (!list1.empty() && !list2.empty()) {
        if (list1.top() == list2.top()) {
            identicalCount++;
            list1.pop();
            list2.pop();
        } else if (list1.top() > list2.top()) {
            list1.pop();
        } else if (list1.top() < list2.top()) {
            list2.pop();
        }
    }
    if (identicalCount > list1Size * lengthRatio * 0.5) {
        return true;
    } else {
        return false;
    }
}

void SeqIterator::getMinHashList(priority_queue <uint64_t> &sortedHashQue, const char *seq) {
    size_t seqLength = strlen(seq);
    size_t kmerLegnth = 24;
    char *kmer = (char *) malloc(sizeof(char) * (kmerLegnth + 1));
    kmer[kmerLegnth] = '\0';
    size_t queLength = 0;
    size_t maxLength = 3000;
    size_t currHash;
    sortedHashQue.push(UINT64_MAX);

    for (size_t i = 0; i + kmerLegnth - 1 < seqLength; i++) {
        strncpy(kmer, seq + i, kmerLegnth);
        currHash = XXH64(kmer, kmerLegnth, 0);
        if (currHash < sortedHashQue.top()) {
            if (queLength < maxLength) {
                sortedHashQue.push(currHash);
                queLength++;
            } else {
                sortedHashQue.pop();
                sortedHashQue.push(currHash);
            }
        }
    }
    free(kmer);
}


void SeqIterator::maskLowComplexityRegions(const unsigned char *seq, unsigned char *maskedSeq, ProbabilityMatrix & probMat,
                                           float maskProb, const BaseMatrix * subMat) {
    unsigned int seqLen = 0;
    while (seq[seqLen] != '\0') {
        maskedSeq[seqLen] = (char) subMat->aa2num[static_cast<int>(seq[seqLen])];
        seqLen++;
    }
    tantan::maskSequences(maskedSeq,
                          maskedSeq + seqLen,
                          50 /*options.maxCycleLength*/,
                          probMat.probMatrixPointers,
                          0.005 /*options.repeatProb*/,
                          0.05 /*options.repeatEndProb*/,
                          0.9 /*options.repeatOffsetProbDecay*/,
                          0, 0,
                          maskProb /*options.minMaskProb*/,
                          probMat.hardMaskTable);
    for (unsigned int pos = 0; pos < seqLen; pos++) {
        char nt = seq[pos];
        maskedSeq[pos] = (maskedSeq[pos] == probMat.hardMaskTable[0]) ? 'N' : nt;
    }
}




void SeqIterator::devideToCdsAndNonCds(const char *maskedSeq,
                                       size_t seqLen,
                                       const vector<CDSinfo> &cdsInfo, 
                                       vector<string> &cds,
                                       vector<string> &nonCds) {
    string tmpMasked;
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        size_t locNum = cdsInfo[i].loc.size();
        int currStartCodonPos = 0;
        for (size_t j = 0; j < locNum; j++) {
            // Extend 21 bases to both sides for k-mer from CDS boudaries
            size_t begin = cdsInfo[i].loc[j].first - 1;
            size_t end = cdsInfo[i].loc[j].second - 1;
            if (j == 0) {
                int k = 0;
                while (k < (this->kmerLen - 1) && begin >= 3) {
                    begin -= 3;
                    currStartCodonPos += 3;
                    k++;
                }
            }
            if (j == locNum - 1) {
                int k = 0;
                while (k < (this->kmerLen - 1) && end + 3 < seqLen) {
                    end += 3;
                    k++;
                }
            }    
            // cout << seqLen << " begin: " << begin << " end: " << end << endl;
            tmpMasked += string(maskedSeq + begin, end - begin + 1);    
        }
        // Reverse complement if needed
        if (cdsInfo[i].isComplement) {
            cds.emplace_back(reverseComplement(tmpMasked));
        } else {
            cds.emplace_back(tmpMasked);
        }
        tmpMasked.clear();
    }

    // Get non-CDS that are not in the CDS list
    bool * checker = new bool[seqLen];
    memset(checker, true, seqLen * sizeof(bool));
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        for (size_t j = 0; j < cdsInfo[i].loc.size(); j++) {
            for (size_t k = cdsInfo[i].loc[j].first - 1; k < cdsInfo[i].loc[j].second; k++) {
                checker[k] = false;
            }
        }
    }

    size_t len;
    size_t i = 0;
    while (i < seqLen) {
        len = 0;
        while (i < seqLen && checker[i]) {
            i++;
            len++;
        }
        if (len > 32) {
            nonCds.emplace_back(string(maskedSeq + i - len, len));
        }
        i ++;
    }
}