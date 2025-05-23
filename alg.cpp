
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include <math>
#include "read_fasta.cpp"

/*
Matrix: m rows, n columns
*/

typedef std::string str;
typedef std::pair<int, int> pair;
typedef std::vector<pair> alignment; // C = ((i_1, j_1), (i_2, j_2), ..., (i_n, j_n)) the indexes of the match (0 if it's a gap)
int INF = std::numeric_limits<int>::min();
int SUP = std::numeric_limits<int>::max();

const int gap_creation_penalty = 2; // h
const int gap_penalty = 1;  // g
const int match_score = 2;

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};

class cell{
public:
    int value;
    struct cell *next, *prev;
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
    if (P.s > 0 and a_type(C.front()) != P.s){return SUP;} // infinity
    return 0;
}

/*
Compute Score Extended Sequence Alignment End
*/
int esce(alignment C, extended_P P){
    if ((P.e == -2 or P.e == -1) and a_type(C.back()) == P.e ){return gap_creation_penalty;};
    if (P.e > 0 and a_type(C.back()) != P.e){return SUP;} // infinity
    return 0;
}

/*
Compute Score Extended Sequence Alignment
*/
int esc(alignment C, extended_P P){
    return escs(C, P)+ esce(C,P) + sc(C,P);
}

/*

TO DO: Add computation of prev and next

*/
void compute_T(int j_first, int j_last, int i_first, int i_last,
               const str& A, const str& B, int s, int e,
               std::function<int(char, char)> f, cell* result1, cell* result2, cell* result3) {

    int l = j_last - j_first + 1;

    std::vector<int> T1_prev(l + 1, INF);
    std::vector<int> T2_prev(l + 1, INF);
    std::vector<int> T3_prev(l + 1, INF);

    std::vector<int> T1_curr(l + 1, INF);
    std::vector<int> T2_curr(l + 1, INF);
    std::vector<int> T3_curr(l + 1, INF);

    // Initialization for row i = i_first - 1 (virtual row above first)
    T1_prev[0] = (s <= 1 ? 0 : INF);

    for (int j = 1; j <= l; j++) {
        if (s == -2){T2_prev[j] = -g * j;}
        else if (s == 2 || s < 0){T2_prev[j] = -h - g * j;}
        if (s == -3){T3_prev[0] = -g;}
        else if (s == 3 || s < 0){T3_prev[0] = -h - g;}
    }

    // Fill DP rows from i_first to i_last
    for (int i = i_first; i <= i_last; i++) {
        // j = 0 column
        if (s == -3){T3_curr[0] = -g * (i - i_first + 1);}
        else if (s == 3 || s < 0){T3_curr[0] = -h - g * (i - i_first + 1);}

        // Row computation: left to right
        for (int j = 1; j <= l; j++) {
            char ai = A[i - 1];       // zero-indexed
            char bj = B[j_first + j - 2]; // adjust index to B

            // T1: Match
            T1_curr[j] = f(ai, bj) + max({
                T1_prev[j - 1],
                T2_prev[j - 1],
                T3_prev[j - 1]
            });

            // T2: Gap in A
            T2_curr[j] = max({
                T1_curr[j - 1] - (g + h),  // Gap opening
                T2_curr[j - 1] - g,        // Gap extension
                T3_curr[j - 1] - (g + h)   // Switch gap type
            });

            // T3: Gap in B
            T3_curr[j] = max({
                T1_prev[j] - (g + h),
                T2_prev[j] - (g + h),
                T3_prev[j] - g
            });
        }

        // Move current row to previous row
        T1_prev = T1_curr;
        T2_prev = T2_curr;
        T3_prev = T3_curr;
    }

    for (int j = j_first; j <= j_last; j++){
        result1[j] = T1_curr[j];
        result2[j] = T2_curr[j];
        result3[j] = T3_curr[j];
    }
    return;
}

void compute_TRev(int j_first, int j_last, int i_first, int i_last,
               const str& A, const str& B, int s, int e,
               std::function<int(char, char)> f, cell* result1, cell* result2, cell* result3) {

    int l = j_last - j_first + 1;

    // T1, T2, T3 represent the next row
    std::vector<int> T1_next(l + 1, INT_MIN);
    std::vector<int> T2_next(l + 1, INT_MIN);
    std::vector<int> T3_next(l + 1, INT_MIN);

    // T1, T2, T3 for the current row
    std::vector<int> T1_curr(l + 1, INT_MIN);
    std::vector<int> T2_curr(l + 1, INT_MIN);
    std::vector<int> T3_curr(l + 1, INT_MIN);

    // Initialization for row i = i_last + 1 (virtual row below last)
    T1_next[l] = (e <= 1 ? 0 : INT_MIN);  // Based on end_type

    for (int j = 0; j < l; j++) {
        if (e == -2){T2_next[j] = -g * (l - j);}
        else if (e == 2 || e < 0){T2_next[j] = -h - g * (l - j);

        if (e == -3){T3_next[l] = -g;}
        else if (e == 3 || e < 0){T3_next[l] = -h - g;}
    }

    // Fill DP rows from i_last to i_first (backwards)
    for (int i = i_last; i >= i_first; i--) {
        // j = l+1 column (after end)
        if (e == -3){T3_curr[l] = -g * (i_last - i + 1);}
        else if (e == 3 || e < 0){T3_curr[l] = -h - g * (i_last - i + 1);}

        // Row computation: right to left
        for (int j = l - 1; j >= 0; j--) {
            char ai = A[i];                       // going forward in A
            char bj = B[j_first + j];             // j is already adjusted

            // T1^R: match/mismatch
            T1_curr[j] = f(ai, bj) + max({
                T1_next[j + 1],
                T2_next[j + 1],
                T3_next[j + 1]
            });

            // T2^R: gap in A
            T2_curr[j] = max({
                T1_curr[j + 1] - (g + h),  // open gap
                T2_curr[j + 1] - g,        // extend gap
                T3_curr[j + 1] - (g + h)   // switch gap
            });

            // T3^R: gap in B
            T3_curr[j] = max({
                T1_next[j] - (g + h),
                T2_next[j] - (g + h),
                T3_next[j] - g
            });
        }

        // Move current row to next row
        T1_next = T1_curr;
        T2_next = T2_curr;
        T3_next = T3_curr;
    }

    for (int j = j_first; j <= j_last; j++){
        result1[j] = T1_curr[j];
        result2[j] = T2_curr[j];
        result3[j] = T3_curr[j];
    }
    return;
}


int decompose(int row_u, int row_d, int col_l, int col_r, int p, int n, int m, int s, int e){
    if ((col_r - col_l < p) or (row_d - row_u) < n/p){STOP?}
    int row_m = ceil((row_u + row_d)/2)*(m/p);t

    int l;
    std::vector<cell> T1(l+1), *T2(l+1), *T3(l+1), *TR1(l+1), *TR2(l+1), *TR3(l+1);
    // compute T and T rev

    // compute opt()
    pair best(INF, 0);

    for (j = col_l; j<col_r;j++){
        int T_max = max(T1[j], T2[j], T3[j]);
        int TR_max = max(TR1[j], TR2[j], TR3[j]);
        int opt = max(T_max + TR_max, T2[j] + TR2[j] + h, T3[j] + TR3[j] + h);
        if (best.first < opt){best = pair(opt, j);}
    }

    //found two special columns and special rows with next prev

    // Call Recursively with 2 threads

    // Solve Partial Problem for this square  (backtrace and score)

    // Join solutions ( add up scores and concatenate strings)

    // return (score, string)
}


// ------------------- COMPUTE ALIGNEMENT SCORE -------------------
int solve(str A, str B, int p){
    int n = A.length();
    int m = B.length();

    decompose(0, m, p, n, m, -1, -1);

}