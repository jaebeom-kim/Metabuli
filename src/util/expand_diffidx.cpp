#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"

using namespace std;

int expand_diffidx(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    // setExpandDiffIdxDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string diffIdxFileName = par.filenames[0];
    string expandedDiffIdxFileName = diffIdxFileName + ".expanded";
    size_t bufferSize = 1000'000'000;

    // Diff idx
    FILE * diffIdxFp = fopen(diffIdxFileName.c_str(), "rb");
    uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (bufferSize + 1)); // size = 32 Mb
    size_t diffIdxBufferIdx = 0;
    size_t diffIdxPos = 0;
    size_t numOfDiffIdx = FileUtil::getFileSize(diffIdxFileName) / sizeof(uint16_t);
    // size_t diffIdxCnt = 0;


    // Expanded idx
    FILE * expandedDiffIdxFp = fopen(expandedDiffIdxFileName.c_str(), "wb");
    uint64_t * expandedIdxBuffer = (uint64_t *) malloc(sizeof(uint64_t) * (bufferSize + 1)); // size = 80 Mb
    size_t expandedIdxBufferIdx = 0;
    fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
    loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, bufferSize);

    bool complete = false;
    uint64_t currentTargetKmer = 0;
    while (!complete) {
        // Load diff idx buffer
        // fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
        size_t loadedBytes = loadBuffer(diffIdxFp,
                                    diffIdxBuffer,
                                diffIdxBufferIdx,
                                    bufferSize);

        // loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                                                   
        if (loadedBytes != bufferSize * sizeof(uint16_t)) {
            complete = true;
        }
        // Expand diff idx in buffer
        while (diffIdxPos != numOfDiffIdx) {
            if (bufferSize < diffIdxBufferIdx + 7) {
                continue;
            }
            currentTargetKmer = KmerMatcher::getNextTargetKmer(currentTargetKmer,
                                                                              diffIdxBuffer,
                                                                               diffIdxBufferIdx,
                                                                                diffIdxPos);
            // diffIdxCnt += 1;
            if (expandedIdxBufferIdx == bufferSize) {
                fwrite(expandedIdxBuffer, sizeof(uint64_t), bufferSize, expandedDiffIdxFp);
                expandedIdxBufferIdx = 0;
            }
            expandedIdxBuffer[expandedIdxBufferIdx++] = currentTargetKmer;
        }
        if (diffIdxPos == numOfDiffIdx) {
            break;
        }
    }
    fwrite(expandedIdxBuffer, sizeof(uint64_t), expandedIdxBufferIdx, expandedDiffIdxFp);
                    
            

    
    return 0;
}

