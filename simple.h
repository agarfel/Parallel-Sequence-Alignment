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
class Entry {
public: 

    int value;
    alignment al;

    Entry() : value(0), al() {}

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