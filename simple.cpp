#include "simple.h"
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include "read_fasta.cpp"
#include <condition_variable>
#include <cmath>
#include <thread>
#include <cstring>
#include <chrono>

int f(char a, char b){
    if(a == b){return 2;} 
    return 0;
}

int from(int best, int t1, int t2, int t3){
    if(t1 == best){return 1;}
    if(t2 == best){return 2;}
    return 3;
}

pair find_last(alignment al){
    pair mypair;
    for(int i = 0; i < al.data.size(); i++){
        if(al.data[i].first != -1){
            mypair.first = al.data[i].first;
            break;
        }

    }
    for(int j = 0; j < al.data.size(); j++){
        if(al.data[j].second != -1){
            mypair.second = al.data[j].second;
            break;
        }
    }

    return mypair;
}

void output_alignement(alignment alignment, str A, str B){
    std::string al_A = "";
    std::string al_B = "";
    int len = alignment.data.size();

    for(int k = 0; k < len; k++){
        int i = alignment.data[k].first;
        if(i == -1){
            al_A += "-";
        }
        else{
            al_A += A[i];
        }
        int j = alignment.data[k].second;
        if(j == -1){
            al_B += "-";
        }
        else{
            al_B += B[j];
        }
    }
    std::cout << "Length Sequence A : " << al_A.size() << std::endl;
    std::cout << "Length Sequence B : " << al_A.size() << std::endl;

    std::cout << "Sequence A : " << al_A << std::endl;
    std::cout << "Sequence B : " << al_B << std::endl;

}


/*
Function called by each thread.
    A_ptr : pointer to A (in global memory)
    B_ptr : pointer to B (in global memory)
    p : number of threads
    id : number of this thread
    n : length of A
    m : lenght B
    Working : stage of each processor
    working_mutexes
    update : conditional variable (updates to working)
    sharingT : to share values of between processors when computing T
    result : store final result
*/
void solve_subproblem_parallel(str& A, str& B_full, int p, int id, std::vector<int>& working,
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<Entry>>* sharingT,
    result& FinalResult){
    int size = B_full.size();
    int num_blocks = p;
    int block_size = size/p;
    int start = id*block_size;
    int B_len;
    if (id == p-1){
        B_len = start + block_size + (size % num_blocks);
    } else {
        B_len = start + block_size;
    }

    std::vector<char> B(B_full.begin()+start, B_full.begin() + B_len); // Copy B part

    int m = B.size();
    int n = A.size();

    bool ishead = (id == 0);
    std::vector<Entry> row1_T1(m + 1);
    std::vector<Entry> row1_T2(m + 1);
    std::vector<Entry> row1_T3(m + 1);

    std::vector<Entry> row2_T1(m + 1);
    std::vector<Entry> row2_T2(m + 1);
    std::vector<Entry> row2_T3(m + 1);

    // initialise 1st rows
    for(int j = 0; j < m +1; j++){
        row1_T1[j].value = -1;
        row1_T3[j].value = -(h + g * j);
    }

    // compute all rows 
    for (int i = 1; i < n+1; i++){ 
        // Get T[i][0]:
        if (ishead){
            row1_T1[0].value = -1;
            row2_T1[0].value = -1; 
            if(i == 1){row1_T1[0].value = 0;}
            
            row1_T2[0].value = -(h + g * (i-1));
            row2_T2[0].value = -(h + g * i);

            row1_T3[0].value = -1;
            row2_T3[0].value = -1;

        } else {
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id]);

                while (working[id] != 2){

                    update.wait(lock);
                }
                row1_T1[0] = (*sharingT)[0][id];
                row1_T2[0] = (*sharingT)[1][id];
                row1_T3[0] = (*sharingT)[2][id];

                row2_T1[0] = (*sharingT)[3][id];
                row2_T2[0] = (*sharingT)[4][id];
                row2_T3[0] = (*sharingT)[5][id];

                working[id] = 1;
            }
            update.notify_all();
        }
        // Compute T_k[i]:
        int t1, t2, t3;
        for(int j = 1; j<m+1; j++){
            // T1
            t1 = row1_T1[j-1].value;
            t2 = row1_T2[j-1].value;
            t3 = row1_T3[j-1].value;
            if ((t1 >= t2) && (t1 >= t3)){ // from T1
                row2_T1[j].value = f(A[i-1], B[j-1]) + t1;
                row2_T1[j].al = row1_T1[j-1].al;
                row2_T1[j].al.data.push_back({i-1, j-1 + start});
            } else  if (t2 >= t3){ // from T2
                row2_T1[j].value =  f(A[i-1], B[j-1]) + t2;
                row2_T1[j].al = row1_T2[j-1].al;
                row2_T1[j].al.data.push_back({-1, j-1 + start});
            } else { // from T3
                row2_T1[j].value =  f(A[i-1], B[j-1]) + t3;
                row2_T1[j].al = row1_T3[j-1].al;
                row2_T1[j].al.data.push_back({i-1, -1});
            }
            // T3
            t1 = row1_T1[j].value - (g+h);
            t2 = row1_T2[j].value - (g+h);
            t3 = row1_T3[j].value - g;
            if ((t1 >= t2) && (t1 >= t3)){ // from T1
                row2_T3[j].value = t1;
                row2_T3[j].al = row1_T1[j].al;
                row2_T3[j].al.data.push_back({i-1, j-1 + start});

            } else  if (t2 >= t3){ // from T2
                row2_T3[j].value = t2;
                row2_T3[j].al = row1_T2[j].al;
                row2_T3[j].al.data.push_back({-1, j-1 + start});

            } else { // from T3
                row2_T3[j].value = t3;
                row2_T3[j].al = row1_T3[j].al;
                row2_T3[j].al.data.push_back({i-1, -1});
            }
            // T2
            t1 = row2_T1[j-1].value - (g+h);
            t2 = row2_T2[j-1].value - g;
            t3 = row2_T3[j-1].value - (g+h);
            if ((t1 >= t2) && (t1 >= t3)){ // from T1
                row2_T2[j].value = t1;
                row2_T2[j].al = row2_T1[j-1].al;
                row2_T2[j].al.data.push_back({i-1, j-1 + start});

            } else if (t2 >= t3){ // from T2
                row2_T2[j].value = t2;
                row2_T2[j].al = row2_T2[j-1].al;
                row2_T2[j].al.data.push_back({-1, j-1 + start});

            } else { // from T3
                row2_T2[j].value = t3;
                row2_T2[j].al = row2_T3[j-1].al;
                row2_T2[j].al.data.push_back({i-1, -1});
            }
        }
        // Share last value:
        if (id != p-1){ // Not last column
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id+1]);

                while (working[id+1] != 1){
                    update.wait(lock);
                }

                (*sharingT)[0][id+1] = row1_T1[m];
                (*sharingT)[1][id+1] = row1_T2[m];
                (*sharingT)[2][id+1] = row1_T3[m];

                (*sharingT)[3][id+1] = row2_T1[m];
                (*sharingT)[4][id+1] = row2_T2[m];
                (*sharingT)[5][id+1] = row2_T3[m];

                working[id+1] = 2;
            }
            update.notify_all();

        }

        std::swap(row1_T1, row2_T1);
        std::swap(row1_T2, row2_T2);
        std::swap(row1_T3, row2_T3);
    }

    // after having computed each row 

    if(id == p-1){ // last processor

        int t1(row1_T1[m].value), t2(row1_T2[m].value), t3(row1_T3[m].value);

        if((t1 >= t2) && (t1 >= t3)){ // best in T1
            FinalResult = result(t1, row1_T1[m].al);
        }
        else if(t2 >= t3){ // best in T2
            FinalResult = result(t2, row1_T2[m].al);
        }
        else { // best in T3
            FinalResult = result(t3, row1_T3[m].al);
        }
        
    }
}

result run(int p, str A, str B){
    std::vector<std::thread> workers(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1);
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<Entry>>* sharingT = new std::vector<std::vector<Entry>>(6, std::vector<Entry>(p));

    result FinalResult;
    for (int i = 0; i<p-1; i++){
        workers[i] = std::thread(solve_subproblem_parallel, std::ref(A), std::ref(B), p,i, std::ref(Working),
        std::ref(working_mutexes), std::ref(update), std::ref(sharingT), std::ref(FinalResult));
    }
    solve_subproblem_parallel(A, B, p, p-1, Working, working_mutexes, update, sharingT, FinalResult);

    for (int i = 0; i < p-1; i++){
        workers[i].join();
    }
   
    pair pair = find_last(FinalResult.al);
    int i(pair.first), j(pair.second);
    std::reverse(FinalResult.al.data.begin(), FinalResult.al.data.end());
    while(i > 0){
        FinalResult.al.data.push_back({i-1, -1});
        i--;
    }
    while(j > 0){
        FinalResult.al.data.push_back({-1, j-1});
        j--;
    }
    std::reverse(FinalResult.al.data.begin(), FinalResult.al.data.end());
    
    return FinalResult;

}


int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <sequence1.fasta> <sequence2.fasta> <p>\n";
        return 1;
    }
    std::string folder = "sequences/";

    std::string sequence1_filename = argv[1];
    std::string sequence2_filename = argv[2];
    int p = std::stoi(argv[3]);

    str A = readFastaSequence(folder + sequence1_filename);
    str B = readFastaSequence(folder+sequence2_filename);

    result Result = run(p, A, B);
    std::cout << "Score: " << Result.score << std::endl;  
    output_alignement(Result.al, A, B);

    // std::cout << "concurrent computation of DP with " << p << " processors" << std::endl;
    // std::cout << "Score of : " << Result.score << std::endl;
    //output_alignement(Result.al, A, B);

    // Uncomment to Test:
    // std::ofstream csvFile("timings.csv");
    // csvFile << "threads,time_microseconds\n";
    
    // for (int p = 1; p <= 32; ++p) {
    //     auto start = std::chrono::high_resolution_clock::now();
    //     
    //     auto end = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //     // Log result
    //     std::cout << "p=" << p << ", time=" << duration.count() << " µs, score=" << Result.score << std::endl;
    //     csvFile << p << "," << duration.count() << "\n";
    // }

    // csvFile.close();
    // std::cout << "Timing results written to timings.csv\n";
    return 0;
}

