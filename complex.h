#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include <condition_variable>
#include <cmath>

typedef std::vector<char> str;
typedef std::pair<int, int> pair;

int INF = -100000;

const int match_score = 2;
const int h = 2; // Gap creation penalty
const int g = 1; // Gap penalty

struct extended_P {
    str A, B;   // Sequences A and B
    int s, e;   // Start and End type for (sub)problem
};

class cell{
public:
    cell() : value(INF) {}
    cell(int v) : value(v) {}
    int value, origin_type, origin_row, r_type, r_row;
};
/*
When aligning substrings of A0 and B0 ,
we still use original positions of characters
in A0 and B0 to model the alignment, as opposed to
using their positions relative to the substrings*/

class alignment {
public:
    std::vector<std::pair<int, int>> data;  // (i,j) i index in A, j index in B (start with i=0 ; -1 ==> Gap)

    alignment(){}

    alignment& operator+=(const alignment& other) {
        data.insert(data.end(), other.data.begin(), other.data.end());
        return *this;
    }

    int type(int k){
        if (k >= data.size()){return -1;}
        if (data[k].first == -1 ){           // Gap in A
            return 2;
        } else if(data[k].second != -1){    // Gap in B
            return 3;
        } else {
            return 1;
        }
    }
};

class result{
public:
    int score;
    alignment al;

    result(){
        score = 0;
        al = {};
    }

    result(int s, alignment a){
        score = s;
        al = a;
    }
    
};


class Info{
public:
    Info(int f, int l, int t, int b, int st, int et, int p=0) {
        first = f;
        last = l;
        top_row = t;
        bottom_row = b;
        s = st;
        e = et;
        plus = p;
    }
    int first, last, top_row, bottom_row, s, e, plus;
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
            data[i] = new int[m]();
        }
    }

~Matrix(){
    for (int i = 0; i < n; i++) {
        delete[] data[i];  // Delete each row
    }
    delete[] data;         // Delete the array of pointers
}


    int* operator[](int i){
        return data[i];
    }
};  


/*
Working values:
    1. Computing T: Read from global
    2. Computing T: Wrote to global
    4. Done computing TRev / Opt
*/