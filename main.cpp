
#include <algorithm>
#include <initializer_list>

#include "read_fasta.cpp"

//---------- Preliminary -------------------------

typedef std::string str;

int match = 1;
int mismatch = -1;
int gap = -2;

class CoupleSequences {
public:

    str seq1, seq2;
    int n, m;

    CoupleSequences(){}

    void load_sequences(const str file1, const str file2){

        this->seq1 = readFastaSequence(file1);
        this->seq2 = readFastaSequence(file2);
        n = seq1.length();
        m = seq2.length();
    }

    int score(int i, int j){
        if(i >= n || j >= m){
            std::cout << "out of range" << std::endl;
            return 0;
        }

       return (seq1[i] == seq2[j]) ? match : mismatch;
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
        for(int i = 0; i<m; i++){
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

//---------- constructing the DP matrix ----------

// Needleman-Wunsh (NW)

void construct_NW(CoupleSequences& Seq, Matrix* H){
    // n rows, m columns
    int n = H->n;
    int m = H->m;

    //initialize first row and first column
    for(int j = 0; j < m; j++){
        (*H)[0][j] = j * gap;
    }
    for(int i = 0; i < n; i++){
        (*H)[i][0] = i * gap;
    }

    //recurrence relation
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            (*H)[i][j] = std::max({((*H)[i-1][j-1] + Seq.score(i,j)), (*H)[i][j-1] + gap, (*H)[i-1][j] + gap});
        }
    }

}

// Gotoh

void construct_Gotoh(int** DP){

    

}

//---------- other stuff -------------------------

int main(){

    CoupleSequences Seq;
    Seq.load_sequences("sequences/Q9CD83.fasta", "sequences/Q9CD83.fasta");
    Matrix* H = new Matrix(Seq.n, Seq.m);
    construct_NW(Seq, H);
    return 0;
}