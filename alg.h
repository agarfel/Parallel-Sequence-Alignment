
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include "read_fasta.cpp"

typedef std::string str;
typedef std::vector<std::pair<int, int>> alignment; // C = ((i_1, j_1), (i_2, j_2), ..., (i_n, j_n)) the indexes of the match (0 if it's a gap)

const int gap_creation_penalty = 2;
const int gap_penalty = 1;
const int match_score = 2;

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};


int a_type(std::pair<int,int>);
int sc(alignment, extended_P);
int escs(alignment, extended_P);
int esce(alignment, extended_P);
int esc(alignment, extended_P);