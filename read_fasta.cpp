#include <string>
#include <fstream>
#include <iostream>


std::string readFastaSequence(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return "";
    }
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue; 
        sequence += line;
    }
    return sequence;
}
