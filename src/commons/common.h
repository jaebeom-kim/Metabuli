#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <cstddef>
#include <utility>
#include "LocalParameters.h"
#include "TaxonomyWrapper.h"
#include <iostream>
#include <unordered_set>
#include "FileUtil.h"
#include <cstdint>


#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#define AA(kmer) ((kmer) & ~16777215)

extern const std::string atcg;
extern const std::string iRCT;

struct MappingRes {
    MappingRes(uint32_t queryId, TaxID speicesId, float score) 
        : queryId(queryId), speciesId(speicesId), score(score) {}
    MappingRes() : queryId(0), speciesId(0), score(0.0f) {}
    uint32_t queryId;
    TaxID speciesId;
    float score; 
};

struct Assembly {
    std::string name;
    TaxID taxid;
    TaxID speciesId;
    TaxID genusId;
    TaxID familyId;
    TaxID orderId;
    Assembly(std::string name) : name(name) {}
    Assembly() : name(""), taxid(0), speciesId(0), genusId(0), familyId(0), orderId(0) {}
    

    void print () const {
        std::cout << "Assembly: " << name << ", TaxID: " << taxid
                  << ", SpeciesID: " << speciesId
                  << ", GenusID: " << genusId
                  << ", FamilyID: " << familyId
                  << ", OrderID: " << orderId << std::endl;
    }
};

struct KmerCnt {
    KmerCnt(size_t length, size_t kmerCnt, size_t totalCnt) : length(length), kmerCnt(kmerCnt), totalCnt(totalCnt) {}
    KmerCnt() : length(0), kmerCnt(0), totalCnt(0){}
    size_t length;
    size_t kmerCnt;
    size_t totalCnt;
};

struct CDSinfo{
    uint32_t protId; //4,294,967,295 counted from 0
    int frame;
    bool isComplement;
    std::vector<std::pair<size_t, size_t>> loc;
    CDSinfo() = default;
    CDSinfo(uint32_t protId, int frame) : protId(protId), frame(frame) {}
};


struct SequenceBlock {
    SequenceBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printSequenceBlock() {
        std::cout << strand << " " << start << " " << end << std::endl;
    }

    int start;
    int end;
    int strand; //true for forward
};

struct Classification {
    Classification() : taxId(0), length(0), score(0.0) {}
    Classification(const std::string & name, TaxID taxId, int length, double score)
        : name(std::move(name)), taxId(taxId), length(length), score(score) {}
    Classification(const std::string & name, int length) 
        : name(std::move(name)), taxId(0), length(length), score(0.0) {}
    std::string name;
    TaxID taxId;
    int length;
    double score;
};
struct Query {
    int queryId;
    int classification;
    float score;
    float coverage;
    int hammingDist;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    int kmerCnt2;
    bool isClassified;
    bool newSpecies; // 36 byte

    std::string name;
    std::map<TaxID,int> taxCnt; // 8 byte per element
    std::vector<std::pair<TaxID, float>> species2Score;
    // std::vector<float> pathScores;

    bool operator==(int id) const { return queryId == id;}

    Query(int queryId, int classification, float score, float coverage, int hammingDist, int queryLength,
          int queryLength2, int kmerCnt, int kmerCnt2, bool isClassified, bool newSpecies, std::string name)
            : queryId(queryId), classification(classification), score(score), coverage(coverage),
              hammingDist(hammingDist), queryLength(queryLength), queryLength2(queryLength2), kmerCnt(kmerCnt), kmerCnt2(kmerCnt2),
              isClassified(isClassified), newSpecies(newSpecies), name(std::move(name)) {}

    Query() : queryId(0), classification(0), score(0), coverage(0), hammingDist(0), queryLength(0),
              queryLength2(0), kmerCnt(0), kmerCnt2(0), isClassified(false), newSpecies(false) {}
};

template<typename T>
struct Buffer {
    T *buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;

    explicit Buffer(size_t sizeOfBuffer=100) {
        buffer = (T *) calloc(sizeOfBuffer, sizeof(T));
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    ~Buffer() {
        if (buffer) {
            free(buffer);
        }
    };

    size_t reserveMemory(size_t numOfKmer) {
        return __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
    };

    void reallocateMemory(size_t sizeOfBuffer) {
        if (sizeOfBuffer > bufferSize) {
            buffer = (T *) realloc(buffer, sizeof(T) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
        }
    };

    void init() {
        startIndexOfReserve = 0;
        if (buffer) {
            memset(buffer, 0, sizeof(T) * bufferSize);
        }
    }
};

template<typename T>
struct ReadBuffer {
    FILE * fp;
    T * p;
    size_t size;
    size_t capacity;
    T * start;
    T * end;

    explicit ReadBuffer(std::string file, size_t sizeOfBuffer=100) {
        fp = fopen(file.c_str(), "rb");
        if (!fp) {
            std::cerr << "Error opening file: " << file << std::endl;
            exit(EXIT_FAILURE);
        }
        p = (T *) calloc(sizeOfBuffer, sizeof(T));
        capacity = sizeOfBuffer;
        size = 0;
        start = p;
        end = p;
    };

    ~ReadBuffer() {
        if (fp) {
            fclose(fp);
        }
        if (start) {
            free(start);
        }
    };

    size_t loadBuffer(size_t unused = 0) {
        memmove(start,                   // dest
                start + (size - unused), // src
                unused * sizeof(T));     // bytes
        size_t readCount = fread(start + unused, sizeof(T), capacity - unused, fp);
        size = readCount + unused;
        end = start + readCount + unused;
        p = start;
        return readCount + unused;
    }

    inline T getNext() {
        if (p >= end) {
            if(loadBuffer() == 0) {
               return T(); // Return default value if no more data 
            }
        }
        return *p++;
    }

    size_t loadBufferAt(size_t offset) {
        fseek(fp, offset * sizeof(T), SEEK_SET);
        size_t readCount = fread(start, sizeof(T), capacity, fp);
        size = readCount;
        end = start + readCount;
        p = start;
        return readCount;
    }
};

template<typename T>
struct WriteBuffer {
    FILE * fp;
    size_t capacity;
    size_t size;
    T * p;
    T * start;
    size_t writeCnt;
    // T * end;

    explicit WriteBuffer(std::string file, size_t sizeOfBuffer=100) {
        fp = fopen(file.c_str(), "wb");
        if (!fp) {
            std::cerr << "Error opening file: " << file << std::endl;
            exit(EXIT_FAILURE);
        }
        p = (T *) calloc(sizeOfBuffer, sizeof(T));
        capacity = sizeOfBuffer;
        size = 0;
        writeCnt = 0;
        start = p;
        // end = p;
    };

    ~WriteBuffer() {
        if (fp) {
            fclose(fp);
        }
        if (start) {
            free(start);
        }
    };

    void flush() {
        fwrite(start, sizeof(T), (p - start), fp);
        p = start;
        size = 0;
    }

    void write(T *data, size_t dataNum = 1) {
        if (size + dataNum > capacity) {
            flush();
        }
        memcpy(p, data, dataNum * sizeof(T));
        p += dataNum;
        size += dataNum;
        writeCnt += dataNum;
    }
};

inline bool fileExist(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void process_mem_usage(double& vm_usage, double& resident_set);

TaxonomyWrapper * loadTaxonomy(const std::string & dbDir, const std::string & taxonomyDir = "");

int loadDbParameters(LocalParameters & par, const std::string & dbDir);

int searchAccession2TaxID(const std::string & name, const std::unordered_map<std::string, int> & acc2taxid);

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number, int cnt) {
    bufferIdx = 0;
    fseek(fp, cnt * sizeof(T), SEEK_CUR);
    return fread(buffer, sizeof(T), number, fp);
}

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number) {
    bufferIdx = 0;
    return fread(buffer, sizeof(T), number, fp);
}


// template <typename T>
// size_t loadBuffer2(FILE *fp, T *buffer, size_t number, int cnt) {
//     fseek(fp, cnt * sizeof(T), SEEK_CUR);
//     return fread(buffer, sizeof(T), number, fp);
// }

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number, size_t unused) {
    memmove(buffer,
            buffer + (number - unused),
            unused * sizeof(T));
    size_t readCount = fread(buffer + unused, sizeof(T), number - unused, fp);
    return readCount + unused;
}

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number) {
    return fread(buffer, sizeof(T), number, fp);
}


template <typename T>
inline T getElement(
    size_t bufferSize,
    FILE *kmerInfoFp,
    T *infoBuffer,
    size_t &infoBufferIdx) 
{
    if (unlikely(infoBufferIdx >= bufferSize)) {
        loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
}

void getObservedAccessionList(const std::string & fnaListFileName,
                              std::vector<std::string> & fastaList,
                              std::unordered_map<std::string, TaxID> & acc2taxid);

void fillAcc2TaxIdMap(std::unordered_map<std::string, TaxID> & acc2taxid,
                      const std::string & acc2taxidFileName);
                                
bool haveRedundancyInfo(const std::string & dbDir);

#endif //ADCLASSIFIER2_COMMON_H
