#include "MetamerPattern.h"
#include "json.h"

using json = nlohmann::json;

MultiCodePattern::MultiCodePattern(const std::string & customFile) {
    std::ifstream file(customFile);
    if (!file.is_open()) {
        std::cout << "Error: Could not open custom metamer pattern file: " << customFile << std::endl;
    }

    std::string line;
    std::stringstream jsonbuf;
    bool inJson = false;
    while (std::getline(file, line)) {
        if (line == "===BEGIN_CUSTOM_METAMER===") {
            inJson = true;
            continue;   
        }
        if (line == "===END_CUSTOM_METAMER===") {
            inJson = false;
            break;  
        }
        if (inJson) {
            jsonbuf << line << "\n";
        }
    }

    json j = json::parse(jsonbuf.str());

    std::string name = j["name"];
    kmerLen = j["length"];
    codePattern = j["position_codes"].get<std::vector<int>>();

    if (kmerLen != static_cast<int>(codePattern.size())) {
        std::cout << "Error: Length of code pattern does not match the specified length." << std::endl;
    }


    for (const auto& codeEntry : j["codes"]) {

        const auto& codons_json = codeEntry["codons"];

        std::cout << "Loading custom genetic code: " << codeEntry["id"] << std::endl;

        std::vector<std::pair<std::string, std::vector<std::string>>> translationTable;

        for (const auto& entry : codons_json) {
            std::string aaStr = entry[0].get<std::string>();
            std::vector<std::string> codons = entry[1].get<std::vector<std::string>>();
            translationTable.emplace_back(std::string{aaStr[0]}, codons);
        }

        geneticCodes.push_back(std::make_unique<CustomGeneticCode>(translationTable));
    }

    // Validate total bits do not exceed 64 bits
    int totalBits = 0;
    for (size_t i = 0; i < codePattern.size(); ++i) {
        totalBits += geneticCodes[codePattern[i]]->bitPerAA;
        totalBits += geneticCodes[codePattern[i]]->bitPerCodon;
        if (totalBits > 64) {
            std::cout << "Error: Total bits exceed 64 bits." << std::endl;
        }
    }

    init();
}