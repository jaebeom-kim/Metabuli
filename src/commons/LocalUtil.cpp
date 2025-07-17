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

int LocalUtil::getFirstWhiteSpacePos(const std::string &str) {
    for (size_t i = 0; i < str.size(); ++i) {
        if (isspace(int(str[i]))) {
            return i;
        }
    }
    return str.size();
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