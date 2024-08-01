#include "LocalUtil.h"
#include <cctype>


std::string LocalUtil::getQueryBaseName(const std::string & queryPath) {
    std::vector<std::string> splits = Util::split(queryPath, ".");
    std::string baseName;
    int extentionNum = 1;
    if (Util::endsWith(".gz", queryPath)) {
        extentionNum = 2;
    }
    for (size_t i = 0; i < splits.size() - extentionNum; ++i) {
        if (i == splits.size() - extentionNum - 1) {
            baseName += splits[i];
        } else {
            baseName += splits[i] + ".";
        }
    }
    return baseName;
}


void LocalUtil::splitQueryFile(std::vector<SequenceBlock> & sequences, const std::string &queryPath) {
    KSeqWrapper* kseq = nullptr;
    kseq = KSeqFactory(queryPath.c_str());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry & e = kseq->entry;
        sequences.emplace_back(e.headerOffset - 1,
                               e.sequenceOffset + e.sequence.l,
                               e.sequenceOffset + e.sequence.l - e.headerOffset + 2,
                               e.sequence.l);
    }
    delete kseq;
}

int LocalUtil::getMaxCoveredLength(int queryLength) {
    if (queryLength % 3 == 2) {
        return queryLength - 2; 
    } else if (queryLength % 3 == 1) {
        return queryLength - 4; 
    } else {
        return queryLength - 3; 
    }
}

int LocalUtil::getFirstWhiteSpacePos(const std::string &str) {
    for (size_t i = 0; i < str.size(); ++i) {
        if (isspace(int(str[i]))) {
            return i;
        }
    }
    return str.size();
}

void LocalUtil::loadMappingFile_text(const std::string &fileName, std::unordered_map<int, int> & map) {
    std::ifstream ifs(fileName);
    if (!ifs) {
        std::cerr << "Could not open " << fileName << " for reading." << std::endl;
        return;
    }

    int key;
    int value;

    while (ifs >> key >> value) {
        map[key] = value;
    }

    ifs.close();
}

// std::string LocalUtil::getAccessionFromHeader(const std::string &header) {
//     int pos = getFirstWhiteSpacePos(header);
//     std::string accession = header.substr(0, pos);
//     std::vector<std::string> splits = Util::split(accession, ".");
//     if (splits.size() > 1) {
//         accession = splits[0];
//     }
//     return std::stoi(accession.substr(3));
// }

void LocalUtil::exportVector(const std::vector<std::string>& vec, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "Could not open " << filename << " for writing." << std::endl;
        return;
    }

    // Write the number of strings
    size_t size = vec.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Write each string
    for (const auto& str : vec) {
        size_t length = str.size();
        ofs.write(reinterpret_cast<const char*>(&length), sizeof(length));
        ofs.write(str.data(), length);
    }
}


 void LocalUtil::importVector(const std::string& filename, std::vector<std::string> & vec) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        std::cerr << "Could not open " << filename << " for reading." << std::endl;
        return;
    }
    // Read the number of strings
    size_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);

    // Read each string
    for (size_t i = 0; i < size; ++i) {
        size_t length;
        ifs.read(reinterpret_cast<char*>(&length), sizeof(length));
        vec[i].resize(length);
        ifs.read(&vec[i][0], length);
    }
}

void LocalUtil::exportVector(const std::vector<int>& vec, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) {
        std::cerr << "Could not open " << filename << " for writing." << std::endl;
        return;
    }

    // Write the number of integers
    size_t size = vec.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Write the integers
    ofs.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(int));
}


void LocalUtil::importVector(const std::string& filename, std::vector<int>& vec) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs) {
        std::cerr << "Could not open " << filename << " for reading." << std::endl;
        return;
    }

    // Read the number of integers
    size_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    vec.resize(size);

    // Read the integers
    ifs.read(reinterpret_cast<char*>(vec.data()), size * sizeof(int));
}
