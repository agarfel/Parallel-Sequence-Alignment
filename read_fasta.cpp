#include <string>
#include <fstream>
#include <iostream>

// the files used to store the sequences are of format .fasta
// this function allows us to get a string version of these files 

std::vector<char> readFastaSequence(const std::string& filename) {
    std::ifstream file(filename);
    std::string line, sequence;
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        std::vector<char> s;
        return s;
    }
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '>') continue; 
        sequence += line;
    }
    std::vector<char> s;
    s.assign(sequence.begin(), sequence.end());
    return s;
}
