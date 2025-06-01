#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include "read_fasta.cpp"
#include <condition_variable>
#include <cmath>

typedef std::vector<char> str;
typedef std::vector<std::pair<int, int>> alignment; // C = ((i_1, j_1), (i_2, j_2), ..., (i_n, j_n)) the indexes of the match (0 if it's a gap)
typedef std::pair<int, int> pair;

int INF = std::numeric_limits<int>::min();
int SUP = std::numeric_limits<int>::max();

const int gap_creation_penalty = 2; // h
const int gap_penalty = 1;  // g
const int match_score = 2;

const int h = gap_creation_penalty;
const int g = gap_penalty;

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};

class cell{
public:
    cell(int v, int n, int p){value = v; next = n; prev = p;}
    int value;
    int next, prev;
};

class result{
    int TO_DO;
};

class Info{
public:
    Info(int f, int l, int t, int b, int st, int et) {
        first = f;
        last = l;
        top_row = t;
        bottom_row = b;
        s = st;
        e = et;
    }
    int first, last, top_row, bottom_row, s, e;
};


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

/*
Working values:
    1. Computing T: Read from global
    2. Computing T: Wrote to global
    3. Done computing T
    4. Done computing Opt
    
    X. Solving sub-problem
*/