
#include <algorithm>
#include <initializer_list>
#include <vector>

#include "read_fasta.cpp"

//---------- Preliminary -----------------------------------

typedef std::string str;

int match = 1;
int mismatch = -1;
int gap = -1;
int gap_op = -4;
int eps = -1e9;

class CoupleSequences {
public:

    str seq1, seq2;
    int n, m;

    CoupleSequences(){}

    CoupleSequences(const str file1, const str file2){
        load_sequences(file1, file2);
    }

    //load sequences from file
    void load_sequences(const str file1, const str file2){
        input_sequences(readFastaSequence(file1), readFastaSequence(file2));
    }

    //input sequences as strings directly 
    void input_sequences(str seq1, str seq2){
        this->seq1 = seq1;
        this->seq2 = seq2;
        n = seq1.length();
        m = seq2.length();
    }

    //score function for the DP matrix creation
    int score(int i, int j){
        if(i >= n || j >= m){
            std::cout << "out of range" << std::endl;
            return 0;
        }

       return (seq1[i] == seq2[j]) ? match : mismatch;
    }

    void print(){
        std::cout << "\nObject Couple of Sequences (" << n << ", " << m << ") \n\n" << "seq1 : " << seq1 << "\n\n" << "seq2 : " << seq2 << "\n" << std::endl;
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

//---------- constructing the DP matrix and traceback ------

// Needleman-Wunsh (NW)

int construct_NW(CoupleSequences& Seq, Matrix& H){
    // n rows, m columns
    int n = H.n;
    int m = H.m;

    //initialize first row and first column
    for(int j = 0; j < m; j++){
        H[0][j] = j * gap;
    }
    for(int i = 0; i < n; i++){
        H[i][0] = i * gap;
    }

    //recurrence relation
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            H[i][j] = std::max({(H[i-1][j-1] + Seq.score(i,j)), H[i][j-1] + gap, H[i-1][j] + gap});
        }
    }
    return H[n-1][m-1];
}

CoupleSequences trace_back_NW(CoupleSequences& Seq, Matrix& H){

    str seq1 = "";
    str seq2 = "";
    int i = H.n - 1;
    int j = H.m - 1;

    while(i > 0 && j > 0){
        int score = H[i][j];

        //comes from top -> gap from seq2
        if(score - gap == H[i-1][j]){
            seq1 = Seq.seq1[i] + seq1;
            seq2 = "-" + seq2;
            i--;
        }
        //comes from left -> gap from seq1
        else if (score - gap == H[i][j-1]){
            seq1 = "-" + seq1;
            seq2 = Seq.seq2[j] + seq2;
            j--;
        }
        //comes from diagonal -> match/mismatch
        else{
            seq1 = Seq.seq1[i-1] + seq1;
            seq2 = Seq.seq2[j-1] + seq2;
            i--;
            j--;
        }
    }

    //leftovers from seq1
    while(i >= 0){
        seq1 = Seq.seq1[i] + seq1;
        seq2 = "-" + seq2;
        i--;
    }

    //leftovers from seq2
    while(j >= 0){
        seq1 = "-" + seq1;
        seq2 = Seq.seq2[j] + seq2;
        j--;
    }

    CoupleSequences Alignement;
    Alignement.input_sequences(seq1, seq2);
    return Alignement;
}
// Gotoh

int construct_Gotoh(CoupleSequences& Seq, Matrix& H, Matrix& E, Matrix& F){
    // n rows, m columns
    int n = H.n;
    int m = H.m;

    //initialize first row and first column
    for(int j = 0; j < m; j++){
        H[0][j] = gap_op + (j-1)*gap;
        F[0][j] = gap_op + (j-1)*gap;
        E[0][j] = eps;
    }
    for(int i = 0; i < n; i++){
        H[i][0] = gap_op + (i-1)*gap;
        E[i][0] = gap_op + (i-1)*gap;
        F[i][0] = eps;
    }

    //recurrence relation
    for(int i = 1; i < n; i++){
        for(int j = 1; j < m; j++){
            E[i][j] = std::max({(E[i][j-1] + gap), (H[i][j-1] + gap_op)});
            F[i][j] = std::max({(F[i-1][j] + gap), (H[i-1][j] + gap_op)});
            H[i][j] = std::max({E[i][j], F[i][j], (H[i-1][j-1] + Seq.score(i, j))}); // 0 in the set ?
        }
    }
    return H[n-1][m-1];

}

CoupleSequences trace_back_Gotoh(CoupleSequences& Seq, Matrix& H, Matrix& E, Matrix& F){

    str seq1 = "";
    str seq2 = "";
    int i = H.n - 1;
    int j = H.m - 1;
    int curr = 0;

    while(i > 0 && j > 0){
        int score;
        
        // In H
        if(curr == 0){
            score = H[i][j];
            seq1 = Seq.seq1[i-1] + seq1;
            seq2 = Seq.seq2[j-1] + seq2;

            //comes from E
            if(score == E[i][j]){
                curr = 1;
            }
            //comes from F
            else if (score == F[i][j]){
                curr = 2;
            }
            //comes from H doesnt change curr

            i--;
            j--;
            continue;
        }

        // In E
        if(curr == 1){
            score = E[i][j];
            seq1 = "-" + seq1;
            seq2 = Seq.seq2[j-1] + seq2;

            //comes from H
            if(score - gap_op == H[i][j-1]){
               curr = 0;
            }
            //comes from E doesnt change curr

            j--;
            continue;
        }
        
        // In F
        if(curr == 2){
            score = F[i][j];
            seq1 = Seq.seq1[i-1] + seq1;
            seq2 = "-" + seq2;

            //comes from H
            if (score - gap_op == H[i-1][j]){
                curr = 0;
            }
            //comes from F doesnt change curr

            i--;
            continue;
        }
        
    }

    //leftovers from seq1
    while(i >= 0){
        seq1 = Seq.seq1[i] + seq1;
        seq2 = "-" + seq2;
        i--;
    }

    //leftovers from seq2
    while(j >= 0){
        seq1 = "-" + seq1;
        seq2 = Seq.seq2[j] + seq2;
        j--;
    }

    CoupleSequences Alignement;
    Alignement.input_sequences(seq1, seq2);
    return Alignement;
}



//---------- other stuff -----------------------------------

// if "G" is entered the Gotoh algorithm is used, otherwise by default the NW is used

int main(int argc, char* argv[]){
    
    bool flag = true;
    if(argc >= 2){
        if(str(argv[1]) == "G"){
            flag = false;
        }
    }

    CoupleSequences Seq;
    str file1 = "sequences/insulin_homo.fasta";
    str file2 = "sequences/insulin_bovin.fasta";
    Seq.load_sequences(file1, file2);
    CoupleSequences alignement;
    int score = 0;

    if(flag){
        Matrix* H = new Matrix(Seq.n, Seq.m);
        score = construct_NW(Seq, *H);
        alignement = trace_back_NW(Seq, *H);
    }
    else{   
        Matrix* H = new Matrix(Seq.n, Seq.m);
        Matrix* E = new Matrix(Seq.n, Seq.m);
        Matrix* F = new Matrix(Seq.n, Seq.m);
        score = construct_Gotoh(Seq, *H, *E, *F);
        alignement = trace_back_Gotoh(Seq, *H, *E, *F);
    }

    std::cout << "score : " << score << std::endl;
    alignement.print();

    return 0;
}