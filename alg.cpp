
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>

#include "read_fasta.cpp"

typedef std::string str;
typedef std::pair<int, int> pair;
typedef std::vector<pair> alignment; // C = ((i_1, j_1), (i_2, j_2), ..., (i_n, j_n)) the indexes of the match (0 if it's a gap)

const int gap_creation_penalty = 2;
const int gap_penalty = 1;
const int match_score = 2;

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};

class Matrix {
public:

    int n, m;
    int** data;

    Matrix(){}

    Matrix(int n, int m){
        this->n = n;
        this->m = m;

        data = new int*[n];
        for(int i = 0; i<n; i++){
            data[i] = new int[m];
        }
    }

    ~Matrix(){
        delete[] data;
    }

    int* operator[](int i){
        return data[i];
    }
};  

/*
Get alignement type 
*/
int a_type(pair p){
    int result;
    if(p.first > 0 and p.second > 0){result = 1;} // Match
    else if(p.first == 0){result = 2;} // Gap in A
    else if(p.second == 0){result = 3;} // Gap in B
    else {
        std::cout << "Type not recognised: "  << p.first << " , " << p.second << std::endl;
    }
}


// ------------------- COMPUTE ALIGNEMENT SCORE -------------------
/*
Compute Score Sequence Alignment

Missing gap penalty somewhere
*/
int sc(alignment C, extended_P P){
    int result = 0;
    for (auto p: C){
        if(P.A[p.first] == P.B[p.second]){result += match_score;}
    }
    return result;
}

/*
Compute Score Extended Sequence Alignment Start
*/
int escs(alignment C, extended_P P){
    if ((P.s == -2 or P.s == -1) and a_type(C.front()) == P.s ){return gap_creation_penalty;};
    if (P.s > 0 and a_type(C.front()) != P.s){return std::numeric_limits<int>::max();} // infinity
    return 0;
}

/*
Compute Score Extended Sequence Alignment End
*/
int esce(alignment C, extended_P P){
    if ((P.e == -2 or P.e == -1) and a_type(C.back()) == P.e ){return gap_creation_penalty;};
    if (P.e > 0 and a_type(C.back()) != P.e){return std::numeric_limits<int>::max();} // infinity
    return 0;
}

/*
Compute Score Extended Sequence Alignment
*/
int esc(alignment C, extended_P P){
    return escs(C, P)+ esce(C,P) + sc(C,P);
}