
#include <algorithm>
#include <initializer_list>
#include <vector>

#include "sequential.cpp"

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