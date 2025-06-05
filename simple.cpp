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

int f(char a, char b){
    if(a == b){return 2;} 
    return 0;
};

int from(int best, int t1, int t2, int t3){
    if(t1 == best){return 1;}
    if(t2 == best){return 2;}
    return 3;
}

/* function should be correct, watch out for out of bound indices in the trace back maybe but otherwise it works */
result solve_subproblem(extended_P sub_problem, int col_offset, int row_offset){ // one or two sub_problems

    str A = sub_problem.A;
    str B = sub_problem.B;
    int start_type = sub_problem.s; // -t1 or -t2
    int end_type = sub_problem.e; // -t3 or t2

    result res;

    /* computing the matrices and the score */

    // n rows, m columns
    int n = A.size();
    int m = B.size();
    Matrix T1(n+1,m+1), T2(n+1,m+1), T3(n+1,m+1);
    Matrix BT1(n+1,m+1), BT2(n+1,m+1), BT3(n+1,m+1);

    // Initialize all cells to -1
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            T1[i][j] = -1;
            T2[i][j] = -1;
            T3[i][j] = -1;
        }
    }

    // initialize top-left corner
    int abs_type = std::abs(start_type);
    
    if (abs_type == 1) {
        T1[0][0] = 0;
    } else if (abs_type == 2) {
        T2[0][0] = 0;
    } else if (abs_type == 3) {
        T3[0][0] = 0;
    }

    // initialize first row
    for (int j = 1; j <= m; ++j) {
        if (start_type == 3) {
            T3[0][j] = -(h + g);
        } else if (start_type == -3) {
            T3[0][j] = -g * j;
        } else if (start_type == -1 || start_type == -2) {
            T3[0][j] = -(h + g * j);
        }
        // else T3[0][i] stays -1
    }

    // initialize first column 
    for (int i = 1; i <= n; ++i) {
        if (start_type == 2) {
            T2[i][0] = -(h + g);
        } else if (start_type == -2) {
            T2[i][0] = -g * i;
        } else if (start_type == -1 || start_type == -3) {
            T2[i][0] = -(h + g * i);
        }
        // else T2[j][0] stays -1
    }

    //recurrence relation
    for(int i = 1; i < n+1; i++){
        for(int j = 1; j < m+1; j++){
            int T1_best = std::max({T1[i-1][j-1], T2[i-1][j-1], T3[i-1][j-1]});
            T1[i][j] = f(A[i-1], B[j-1]) + T1_best;
            BT1[i][j] = from(T1_best, T1[i-1][j-1], T2[i-1][j-1], T3[i-1][j-1]);
            int T2_best = std::max({(T1[i][j-1] -(g+h)), (T2[i][j-1] - g), (T3[i][j-1] -(g+h))});
            T2[i][j] = T2_best;
            BT2[i][j] = from(T2_best, (T1[i][j-1] -(g+h)), (T2[i][j-1] - g), (T3[i][j-1] -(g+h)));
            int T3_best = std::max({(T1[i-1][j] -(g+h)), (T2[i-1][j] -(g+h)), (T3[i-1][j] - g)});
            T3[i][j] = T3_best;
            BT3[i][j] = from(T3_best, (T1[i-1][j] -(g+h)), (T2[i-1][j] -(g+h)), (T3[i-1][j] - g));
        }   
    }

    if (end_type > 0) {
        if (end_type == 1){ res.score = T1[n][m];}
        else if (end_type == 2) {res.score = T2[n][m];}
        else {res.score = T3[n][m];}
    } else {
        res.score = std::max({ T1[n][m], T2[n][m], T3[n][m] });
    }

    /* doing the traceback */
    int state = from(res.score, T1[n][m], T2[n][m], T3[n][m]);
    int i(n), j(m);

    while (i > 0 || j > 0) {
        if (state == 1) { // T1
            res.al.data.push_back({i-1 + row_offset, j-1+ col_offset});
            int bt = BT1[i][j];
            i--; j--;
            state = bt;
        } else if (state == 2) { // T2
            res.al.data.push_back({-1, j-1+ col_offset});
            int bt = BT2[i][j];
            j--;
            state = bt;
        } else if (state == 3) { // T3
            res.al.data.push_back({i-1+ row_offset, -1});
            int bt = BT3[i][j];
            i--;
            state = bt;
        }
    }

    std::reverse(res.al.data.begin(), res.al.data.end());
    return res;
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
int solve_subproblem_parallel(str A, str B, int s, int e, int p, int id, std::vector<int>& working,
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<int>>* sharingT,
    result& result){
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
    std::vector<std::vector<int>> T1(num_rows + 1, std::vector<int>(num_cols + 1));
    std::vector<std::vector<int>> T2(num_rows + 1, std::vector<int>(num_cols + 1));
    std::vector<std::vector<int>> T3(num_rows + 1, std::vector<int>(num_cols + 1));

    int m = num_cols;
    for (int j = 0; j < num_cols; j++){          
            if (!((s == 1) || (s == 3))){ T2[0][j] = -h-g*(j + 1);}
            if (s == -2){ T2[0][j] += h;}
    }

    if(ishead){
        if (s == -1){ T1[0][0] = 0; }
        if (s == -2){ T2[0][0] = 0;}
        if (s == -3){ T3[0][0] = 0;}
    }

    for (int i = 1; i < num_rows; i++){ // COMPUTE ALL ROWS
        // Get T[i][0]:
        if (ishead){
            T1[i][0] = INF; // NEEDED?
            T2[i][0] = INF; // NEEDED?
            T3[i][0] = INF;
            if (!((s == 1) || (s == 2))){ T3[i][0] = -h-g*(i + 1);}
            if (s == -3){ T3[i][0] += h;}
        } else {
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id]);

                while (working[id] != 2){

                    update.wait(lock);
                }
                T1[i][0] = (*sharingT)[0][id];
                T2[i][0] = (*sharingT)[1][id];
                T3[i][0] = (*sharingT)[2][id];
            
                working[id] = 1;
            }
            update.notify_all();
        }
        // Compute T[i]:
        for(int j = 1; j<m; j++){
            // T1
            if ((T1[i-1][j-1] >= T2[i-1][j-1]) && (T1[i-1][j-1] >= T3[i-1][j-1])){
                T1[i][j] = f(A[i], B[j]) + T1[i-1][j-1];
            } else  if (T2[i-1][j-1] >= T3[i-1][j-1]){
                T1[i][j] = f(A[i], B[j]) + T2[i-1][j-1];
            } else {
                T1[i][j] = f(A[i], B[j]) + T3[i-1][j-1];
            }
            // T3
            int t1 = T1[i-1][j] - (g+h);
            int t2 = T2[i-1][j] - (g+h);
            int t3 = T3[i-1][j] - g;
            if ((t1 >= t2) && (t1 >= t3)){
                T3[i][j] = t1;

            } else  if (t2 >= t3){
                T3[i][j] = t2;

            } else {
                T3[i][j] = t3; 
            }
            // T2
            t1 = T1[i][j-1] - (g+h);
            t2 = T2[i][j-1] - g;
            t3 = T3[i][j-1] - (g+h);
            if ((t1 >= t2) && (t1 >= t3)){
                T2[i][j] = t1;

            } else if (t2 >= t3){
                T2[i][j] = t2;

            } else {
                T2[i][j] = t3;
            }
        }
        // Share last value:
        if (id != col_k2*2){ // Not last row
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id+2]);

                while (working[id+2] != 1){
                    update.wait(lock);
                }

                (*sharingT)[0][id+2] = T1[i][0];
                (*sharingT)[1][id+2] = T2[i][0];
                (*sharingT)[2][id+2] = T3[i][0];
                working[id+2] = 2;
            }
            update.notify_all();

        }
    }
}

void run(int p, str A, str B){
    std::vector<std::thread> workers(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1); // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<cell>>* sharingT = new std::vector<std::vector<cell>>(3, std::vector<cell>(p));

    for (int i = 0; i<p-1; i++){
        workers[i] = std::thread(solve_subproblem_parallel, A.data(), B.data(), p,i, std::ref(Working),
        std::ref(working_mutexes), std::ref(update), sharingT);
    }
    int TotalScore = solve_subproblem_parallel(A, B, p,p, &Working,
        &working_mutexes, &update, &sharingT);
    alignment FinalAlignment;
    for (int i = 0; i < p; i++){
        workers[i].join();
    }
}