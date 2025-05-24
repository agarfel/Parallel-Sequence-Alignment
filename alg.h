
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include "read_fasta.cpp"
#include <condition_variable>
#include <math>
#include "read_fasta.cpp"

typedef std::string str;
typedef std::vector<std::pair<int, int>> alignment; // C = ((i_1, j_1), (i_2, j_2), ..., (i_n, j_n)) the indexes of the match (0 if it's a gap)
typedef std::pair<int, int> pair;

int INF = std::numeric_limits<int>::min();
int SUP = std::numeric_limits<int>::max();

const int gap_creation_penalty = 2; // h
const int gap_penalty = 1;  // g
const int match_score = 2;

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};

std::vector<int> Working; // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
std::condition_variable update; // notifies if there's a change to Working


class Matrix {
public:

    int n, m;
    struct cell** data;

    Matrix(){}

    Matrix(int n, int m){
        this->n = n;
        this->m = m;

        data = new cell*[n];
        for(int i = 0; i<n; i++){
            data[i] = new cell[m];
        }
    }

    ~Matrix(){
        delete[] data;
    }

    cell* operator[](int i){
        return data[i];
    }
};  

int a_type(std::pair<int,int>);
int sc(alignment, extended_P);
int escs(alignment, extended_P);
int esce(alignment, extended_P);
int esc(alignment, extended_P);

void thread_f(str B, str* A, int first, int last);