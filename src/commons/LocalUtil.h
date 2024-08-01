#ifndef METABULI_LOCALUTIL_H
#define METABULI_LOCALUTIL_H

#include "Util.h"
#include <string>
#include <unordered_map>
#include "common.h"
#include "KSeqWrapper.h"
#include <iostream>
#include <fstream>

class LocalUtil : public Util {
public:
    LocalUtil() = default;

    static std::string getQueryBaseName(const std::string & queryPath);

    template<typename T>
    static T getQueryKmerNumber(T queryLength, int spaceNum);

    template<typename K, typename V>
    static void writeMappingFile(const std::unordered_map<K, V> & map, const std::string & fileName);

    template <typename K, typename V>
    static void loadMappingFile(const std::string& fileName, std::unordered_map<K, V>& map);

    template<typename K, typename V>
    static void writeMappingFile_text(const std::unordered_map<K, V> & map, const std::string & fileName);

    static void loadMappingFile_text(const std::string& fileName, std::unordered_map<int, int>& map);

    static void splitQueryFile(std::vector<SequenceBlock> & seqSegments, const std::string & queryPath);

    static int getMaxCoveredLength(int queryLength) ;

    static int getFirstWhiteSpacePos(const std::string & str);

    static void exportVector(const std::vector<std::string>& vec, const std::string& filename);
    static void exportVector(const std::vector<int>& vec, const std::string& filename); 

    static void importVector(const std::string& filename, std::vector<std::string>& vec);
    static void importVector(const std::string& filename, std::vector<int>& vec);

    static TaxID getTaxIdFromComment(const std::string &header) {
        // comment format: 
        // peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE
        size_t taxIDPos = header.find("TaxID=");
        size_t taxIDEndPos = header.find(" ", taxIDPos);
        return std::stoi(header.substr(taxIDPos + 6, taxIDEndPos - taxIDPos - 6));    
    }


    // static std::string getAccessionFromHeader(const std::string & header);
};


template <typename T>
T LocalUtil::getQueryKmerNumber(T queryLength, int spaceNum) {
    return (getMaxCoveredLength(queryLength) / 3 - kmerLength - spaceNum + 1) * 6;
}

template <typename K, typename V>
void LocalUtil::writeMappingFile(const std::unordered_map<K, V> & map, const std::string & fileName) {
    // Write the mapping files
    std::ofstream ofs(fileName, std::ios::binary);
    if (!ofs) {
        std::cerr << "Could not open " << fileName << " for writing." << std::endl;
        return;
    }

    // Write the size of the map
    size_t size = map.size();
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // Write each key-value pair
    size_t keySize = sizeof(K);
    size_t valueSize = sizeof(V);
    for (const auto& pair : map) {
        ofs.write(reinterpret_cast<const char*>(&pair.first), keySize);
        ofs.write(reinterpret_cast<const char*>(&pair.second), valueSize);
    }
    ofs.close();
}

template <typename K, typename V>
void LocalUtil::loadMappingFile(const std::string &fileName, std::unordered_map<K, V> &map) {
    std::ifstream ifs(fileName, std::ios::binary);
    if (!ifs) {
        std::cerr << "Could not open " << fileName << " for reading." << std::endl;
        return;
    }

    size_t size;
    ifs.read(reinterpret_cast<char*>(&size), sizeof(size));
    std::cout << size << std::endl;

    K key;
    V value;
    size_t keySize = sizeof(K);
    size_t valueSize = sizeof(V);

    for (size_t i = 0; i < size; ++i) {
        ifs.read(reinterpret_cast<char*>(&key), keySize);
        ifs.read(reinterpret_cast<char*>(&value), valueSize);
        map[key] = value;
    }

    ifs.close();
}

template <typename K, typename V>
void LocalUtil::writeMappingFile_text(const std::unordered_map<K, V> & map, const std::string & fileName) {
    // Write the mapping files
    std::ofstream ofs(fileName);
    if (!ofs) {
        std::cerr << "Could not open " << fileName << " for writing." << std::endl;
        return;
    }

    // Write each key-value pair
    for (const auto& pair : map) {
        ofs << pair.first << "\t" << pair.second << std::endl;
    }
    ofs.close();
}




#endif //METABULI_LOCALUTIL_H
