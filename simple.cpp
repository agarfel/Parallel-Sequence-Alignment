#include "alg.h"
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

class Entry {
public: 

    int value;
    alignment al;

    Entry() : value(0), al() {}

};

int f(char a, char b){
    if(a == b){return 2;} 
    return 0;
}

int from(int best, int t1, int t2, int t3){
    if(t1 == best){return 1;}
    if(t2 == best){return 2;}
    return 3;
}



/*
Function called by each thread.
    A_ptr : pointer to B (in global memory)
    B_ptr : pointer to A (in global memory)
    p : number of threads
    id : number of this thread
    n : length of A
    m : lenght B
    Working
    working_mutexes
    update
    sharingT
    sharingOpt
*/
int solve_subproblem_parallel(str A, str B, int p, int id, std::vector<int>& working,
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<Entry>>* sharingT,
    result& FinalResult){
    int m = B.size();
    int num_blocks = p / 2;
    int B_len;
    if ((id == p) || (id == p-1)){
        B_len = (m / num_blocks) + (m % num_blocks);
    } else {
        B_len = m/num_blocks;
    }

    std::vector<char> B(B.begin(), B.begin() + B_len); // Copy B part

    int num_cols = B.size();
    int num_rows = A.size();
    int col_k1 = 0;
    int col_k2 = p;
    int row_k1 = 0;
    int row_k2 = A.size();

    bool ishead = (2*col_k1 == id);
    std::vector<Entry> row1_T1(num_rows + 1);
    std::vector<Entry> row1_T2(num_rows + 1);
    std::vector<Entry> row1_T3(num_rows + 1);

    std::vector<Entry> row2_T1(num_rows + 1);
    std::vector<Entry> row2_T2(num_rows + 1);
    std::vector<Entry> row2_T3(num_rows + 1);

    int m = num_cols;
   
    if(ishead){
        row1_T1[0].value = 0;
    }

    for (int i = 1; i < num_rows; i++){ // COMPUTE ALL ROWS
        // Get T[i][0]:
        if (ishead){
            row1_T1[0].value = INF; // NEEDED?
            row1_T2[0].value = INF; // NEEDED?
           row1_T3[0].value = -h-g*(i + 1);
            row1_T3[0].value += h;
        } else {
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id]);

                while (working[id] != 2){

                    update.wait(lock);
                }
                row2_T1[0] = (*sharingT)[0][id];
                row2_T2[0] = (*sharingT)[1][id];
                row2_T3[0] = (*sharingT)[2][id];
            
                working[id] = 1;
            }
            update.notify_all();
        }
        // Compute T[i]:
        for(int j = 1; j<m; j++){
            // T1
            if ((row1_T1[j-1].value >= row1_T2[j-1].value) && (row1_T1[j-1].value >= row1_T3[j-1].value)){ // from T1
                row2_T1[j].value = f(A[i], B[j]) + row1_T1[j-1].value;
                row2_T1[j].al = row1_T1[j-1].al;
                row2_T1[j].al.data.push_back({i, j});
            } else  if (row1_T2[j-1].value >= row1_T3[j-1].value){ // from T2
                row2_T1[j].value =  f(A[i], B[j]) + row1_T2[j-1].value;
                row2_T1[j].al = row1_T2[j-1].al;
                row2_T1[j].al.data.push_back({-1, j});
            } else { // from T3
                row2_T1[j].value =  f(A[i], B[j]) + row1_T3[j-1].value;
                row2_T1[j].al = row1_T3[j-1].al;
                row2_T1[j].al.data.push_back({i, -1});
            }
            // T3
            int t1 = row1_T1[j].value - (g+h);
            int t2 = row1_T2[j].value - (g+h);
            int t3 = row1_T3[j].value - g;
            if ((t1 >= t2) && (t1 >= t3)){ // from T1
                row2_T3[j].value = t1;
                row2_T3[j].al = row1_T1[j].al;
                row2_T3[j].al.data.push_back({i, j});

            } else  if (t2 >= t3){ // from T2
                row2_T3[j].value = t2;
                row2_T3[j].al = row1_T2[j].al;
                row2_T3[j].al.data.push_back({-1, j});

            } else { // from T3
                row2_T3[j].value = t3;
                row2_T3[j].al = row1_T3[j].al;
                row2_T3[j].al.data.push_back({i, -1});
            }
            // T2
            t1 = row2_T1[j-1].value - (g+h);
            t2 = row2_T2[j-1].value - g;
            t3 = row2_T3[j-1].value - (g+h);
            if ((t1 >= t2) && (t1 >= t3)){ // from T1
                row2_T2[j].value = t1;
                row2_T2[j].al = row2_T1[j-1].al;
                row2_T2[j].al.data.push_back({i, j});

            } else if (t2 >= t3){ // from T2
                row2_T2[j].value = t2;
                row2_T2[j].al = row2_T2[j-1].al;
                row2_T2[j].al.data.push_back({-1, j});

            } else { // from T3
                row2_T2[j].value = t3;
                row2_T2[j].al = row2_T3[j-1].al;
                row2_T2[j].al.data.push_back({i, -1});
            }
        }
        // Share last value:
        if (id != col_k2*2){ // Not last row
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id+2]);

                while (working[id+2] != 1){
                    update.wait(lock);
                }

                (*sharingT)[0][id+2] = row2_T1[0];
                (*sharingT)[1][id+2] = row2_T2[0];
                (*sharingT)[2][id+2] = row2_T3[0];
                working[id+2] = 2;
            }
            update.notify_all();

        }

        std::swap(row1_T1, row2_T1);
        std::swap(row1_T2, row2_T2);
        std::swap(row1_T3, row2_T3);
    }

    // after having computed each row 

    if(id == p-1){ // last processor

        int t1(row1_T1[m-1].value), t2(row1_T2[m-1].value), t3(row1_T3[m-1].value);
        if((t1 >= t2) && (t1 >= t3)){ // best in T1
            FinalResult = result(t1, row1_T1[m-1].al);
        }
        else if(t2 >= t3){ // best in T2
            FinalResult = result(t2, row1_T2[m-1].al);
        }
        else { // best in T3
            FinalResult = result(t3, row1_T3[m-1].al);
        }
        
    }
}

result run(int p, str A, str B){
    std::vector<std::thread> workers(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1); // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<Entry>>* sharingT = new std::vector<std::vector<Entry>>(3, std::vector<Entry>(p));

    result FinalResult;
    for (int i = 0; i<p-1; i++){
        workers[i] = std::thread(solve_subproblem_parallel, A.data(), B.data(), p,i, std::ref(Working),
        std::ref(working_mutexes), std::ref(update), sharingT, std::ref(FinalResult));
    }
    solve_subproblem_parallel(A, B, p, p-1, Working, working_mutexes, update, sharingT, FinalResult);

    for (int i = 0; i < p; i++){
        workers[i].join();
    }

    return FinalResult;

}