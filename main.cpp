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
result solve_subproblem(extended_P sub_problem){ // one or two sub_problems

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
            res.al.data.push_back({i-1, j-1});
            int bt = BT1[i][j];
            i--; j--;
            state = bt;
        } else if (state == 2) { // T2
            res.al.data.push_back({-1, j-1});
            int bt = BT2[i][j];
            j--;
            state = bt;
        } else if (state == 3) { // T3
            res.al.data.push_back({i-1, -1});
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
void thread_f(const char* A_ptr, const char* B_ptr, int p, int id, int n, int m, std::vector<int>& working,
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<cell>>* sharingT,
    std::vector<std::vector<cell>>* sharingTRev, std::vector<cell>* sharingOpt, std::vector<Info>& info,
    result& result){

    int num_blocks = p / 2;
    int B_len;
    if ((id == p) || (id == p-1)){
        B_len = (m / num_blocks) + (m % num_blocks);
    } else {
        B_len = m/num_blocks;
    }

    std::vector<char> B(B_ptr, B_ptr + B_len); // Copy B part

    int num_cols = B.size();
    int num_special_rows = n / p; 
    int col_k1 = info[id].first; // index of first special column (in group)
    int col_k2 = info[id].last; // index of last special column (in group)
    int row_k1 = info[id].top_row; // index of first special row (in group)
    int row_k2 = info[id].bottom_row; // index of last special row (in group)
    int s = info[id].s;
    int e = info[id].e;



    // Decomposing
    while ((col_k1 + 1 != col_k2) && (row_k1 +1 != row_k2)){
        
        int row_midk = ceil((row_k1 + row_k2)/2);   //index of middle row
        bool ishead = (2*col_k1 == id);

        // Copy Section of A we want
        int num_rows = (row_midk - row_k1)*p;
        std::vector<char> A(num_rows);
        std::memcpy(A.data(), A_ptr + row_k1*(n/p), num_rows);

        // Initialize curr and prev
        std::vector<cell> T1_curr(num_cols + 1);
        std::vector<cell> T2_curr(num_cols + 1);
        std::vector<cell> T3_curr(num_cols + 1);

        std::vector<cell> T1_prev(num_cols + 1, INF);
        std::vector<cell> T2_prev(num_cols + 1, INF);
        std::vector<cell> T3_prev(num_cols + 1, INF);

        if(id % 2 == 0){
            int m = num_cols; // amount of columns per thread

            // INITIALIZE T_prev for all threads:
            for (int j = 0; j < num_cols; j++){          
                    if (!((s == 1) || (s == 3))){ T2_prev[j].value = -h-g*(j + 1);}
                    if (s == -2){ T2_prev[j].value += h;}
            }

            if(ishead){
                // [i'-1, j'-1]
                if (s == -1){ T1_prev[0].value = 0; }
                if (s == -2){ T2_prev[0].value = 0;}
                if (s == -3){ T3_prev[0].value = 0;}
            }

            for (int i = 0; i < num_rows; i++){ // COMPUTE ALL ROWS
                // Get T_curr[0]:
                if (ishead){
                    T1_curr[0].value = INF; // NEEDED?
                    T2_curr[0].value = INF; // NEEDED?
                    T3_curr[0].value = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[0].value = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[0].value += h;}
                    T1_curr[0].origin_row = i;
                    T1_curr[0].origin_type = 1;
                    T2_curr[0].origin_row = i;
                    T2_curr[0].origin_type = 2;
                    T3_curr[0].origin_row = i;
                    T3_curr[0].origin_type = 3;
                } else {
                    std::unique_lock<std::mutex> lock(working_mutexes[id]);
                    while (working[id] != 2){
                        update.wait(lock);
                    }
                    T1_curr[0] = (*sharingT)[0][id];
                    T2_curr[0] = (*sharingT)[1][id];
                    T3_curr[0] = (*sharingT)[2][id];
                    {
                        std::lock_guard<std::mutex> lock(working_mutexes[id]);
                        working[id] = 1;
                    }
                    update.notify_all();
                }

                // Compute T_curr:
                for(int j = 1; j<m; j++){
                    // T1
                    if ((T1_prev[j-1].value >= T2_prev[j-1].value) && (T1_prev[j-1].value >= T3_prev[j-1].value)){
                        T1_curr[j].value = f(A[i], B[j]) + T1_prev[j-1].value;
                        T1_curr[j].origin_row = T1_prev[j-1].origin_row;
                        T1_curr[j].origin_type = T1_prev[j-1].origin_type;

                    } else  if (T2_prev[j-1].value >= T3_prev[j-1].value){
                        T1_curr[j].value = f(A[i], B[j]) + T2_prev[j-1].value;
                        T1_curr[j].origin_row = T2_prev[j-1].origin_row;
                        T1_curr[j].origin_type = T2_prev[j-1].origin_type;

                    } else {
                        T1_curr[j].value = f(A[i], B[j]) + T3_prev[j-1].value;
                        T1_curr[j].origin_row = T3_prev[j-1].origin_row;
                        T1_curr[j].origin_type = T3_prev[j-1].origin_type;   
                    }
                    // T3
                    int t1 = T1_prev[j].value - (g+h);
                    int t2 = T2_prev[j].value - (g+h);
                    int t3 = T3_prev[j].value - g;
                    if ((t1 >= t2) && (t1 >= t3)){
                        T3_curr[j].value = t1;
                        T3_curr[j].origin_row = T1_prev[j].origin_row;
                        T3_curr[j].origin_type = T1_prev[j].origin_type;

                    } else  if (t2 >= t3){
                        T3_curr[j].value = t2;
                        T3_curr[j].origin_row = T2_prev[j].origin_row;
                        T3_curr[j].origin_type = T2_prev[j].origin_type;

                    } else {
                        T3_curr[j].value = t3;
                        T3_curr[j].origin_row = T3_prev[j].origin_row;
                        T3_curr[j].origin_type = T3_prev[j].origin_type;   
                    }
                    // T2
                    t1 = T1_curr[j-1].value - (g+h);
                    t2 = T2_curr[j-1].value - g;
                    t3 = T3_curr[j-1].value - (g+h);
                    if ((t1 >= t2) && (t1 >= t3)){
                        T2_curr[j].value = t1;
                        T2_curr[j].origin_row = T1_curr[j-1].origin_row;
                        T2_curr[j].origin_type = T1_curr[j-1].origin_type;

                    } else if (t2 >= t3){
                        T2_curr[j].value = t2;
                        T2_curr[j].origin_row = T2_curr[j-1].origin_row;
                        T2_curr[j].origin_type = T2_curr[j-1].origin_type;

                    } else {
                        T2_curr[j].value = t3;
                        T2_curr[j].origin_row = T3_curr[j-1].origin_row;
                        T2_curr[j].origin_type = T3_curr[j-1].origin_type;   
                    }
                    if (j == 1){
                        T1_curr[1].origin_row = i;
                        T1_curr[1].origin_type = 1;
                        T2_curr[1].origin_row = i;
                        T2_curr[1].origin_type = 2;
                        T3_curr[1].origin_row = i;
                        T3_curr[1].origin_type = 3;
                    }
                }

                // Share last value:
                if (id != col_k2*2){ // Not last row
                    std::unique_lock<std::mutex> lock(working_mutexes[id+2]);
                    while (working[id+2] != 1){
                        update.wait(lock);
                    }
                    (*sharingT)[0][id+2] = T1_curr[0];
                    (*sharingT)[1][id+2] = T2_curr[0];
                    (*sharingT)[2][id+2] = T3_curr[0];
                    {
                        std::lock_guard<std::mutex> lock(working_mutexes[id+2]);
                        working[id+2] = 2;
                    }
                    update.notify_all();
                }

                std::swap(T1_prev, T1_curr);
                std::swap(T2_prev, T2_curr);
                std::swap(T3_prev, T3_curr);

            }

            // Wait until TRev Computed for this set of columns
            std::unique_lock<std::mutex> lock(working_mutexes[id+1]);
            while (working[id+1] != 4){
                update.wait(lock);
            }

            // Copmute opt
            cell max_opt(INF);
            int start = p * id * 2;
            int end = start + num_cols;
            std::vector<cell> T1_Rev(sharingTRev->at(0).begin() + start, sharingTRev->at(0).begin() + end);
            std::vector<cell> T2_Rev(sharingTRev->at(1).begin() + start, sharingTRev->at(1).begin() + end);
            std::vector<cell> T3_Rev(sharingTRev->at(2).begin() + start, sharingTRev->at(2).begin() + end);

            for (int j = 0; j < num_cols; j++){
                int tmp;
                int o_row, r_row, o_type, r_type, TRevmax, Tmax;
                if ((T1_Rev[j].value >= T2_Rev[j].value) && (T1_Rev[j].value >= T3_Rev[j].value)){
                    r_row = T1_Rev[j].r_row;
                    r_type = T1_Rev[j].r_type;
                    TRevmax = T1_Rev[j].value;
                } else if (T2_Rev[j].value >= T3_Rev[j].value){
                    r_row = T2_Rev[j].r_row;
                    r_type = T2_Rev[j].r_type;
                    TRevmax = T2_Rev[j].value;
                } else {
                    r_row = T3_Rev[j].r_row;
                    r_type = T3_Rev[j].r_type;
                    TRevmax = T3_Rev[j].value;
                }
                if ((T1_curr[j].value >= T2_curr[j].value) && (T1_curr[j].value >= T2_curr[j].value)){
                    o_row = T1_curr[j].r_row;
                    o_type = T1_curr[j].r_type;
                    TRevmax = T1_curr[j].value;
                } else if (T2_curr[j].value >= T3_curr[j].value){
                    o_row = T2_curr[j].r_row;
                    o_type = T2_curr[j].r_type;
                    TRevmax = T2_curr[j].value;
                } else {
                    o_row = T3_curr[j].r_row;
                    o_type = T3_curr[j].r_type;
                    TRevmax = T3_curr[j].value;
                }
                tmp = TRevmax + Tmax;
                if (tmp > max_opt.value){
                    max_opt.value = tmp;
                    max_opt.origin_type = o_type;
                    max_opt.r_type = r_type;
                    max_opt.origin_row = o_row;
                    max_opt.r_row = r_row;
                }
            }
            // Share opt
            (*sharingOpt)[id/2] = max_opt;
            {
                std::lock_guard<std::mutex> lock(working_mutexes[id]);
                working[id] = 4;
            }
            update.notify_all();


        } else {
            // Compute TRev
            int m = num_cols; // amount of columns per thread

            // INITIALIZE T_prev for all threads:
            for (int j = num_cols -1; j >= 0; j--){          
                    if (!((e == 1) || (e == 3))){ T2_prev[j].value = -h-g*(j + 1);}
                    if (e == -2){T2_prev[j].value += h;}
            }

            if(id == col_k2+1){ // Is "sub-head"
                // [i'-1, j'-1]
                if (s == -1){ T1_prev[0].value = 0; }
                if (s == -2){ T2_prev[0].value = 0;}
                if (s == -3){ T3_prev[0].value = 0;}
                
            }

            for (int i = num_rows -1; i >= 0; i--){ // COMPUTE ALL ROWS
                // Get T_curr[0]:
                if (id == p){
                    T1_curr[num_cols -1].value = INF; // NEEDED?
                    T2_curr[num_cols -1].value = INF; // NEEDED?
                    T3_curr[num_cols -1].value = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[num_cols -1].value = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[num_cols -1].value += h;}
                    T1_curr[num_cols -1].r_row = i;
                    T1_curr[num_cols -1].r_type = 1;
                    T2_curr[num_cols -1].r_row = i;
                    T2_curr[num_cols -1].r_type = 2;
                    T3_curr[num_cols -1].r_row = i;
                    T3_curr[num_cols -1].r_type = 3;
                } else {
                    std::unique_lock<std::mutex> lock(working_mutexes[id]);
                    while (working[id] != 2){
                        update.wait(lock);
                    }
                    T1_curr[0] = (*sharingT)[0][id];
                    T2_curr[0] = (*sharingT)[1][id];
                    T3_curr[0] = (*sharingT)[2][id];
                    {
                        std::lock_guard<std::mutex> lock(working_mutexes[id]);
                        working[id] = 1;
                    }
                    update.notify_all();
                }

                // Compute T_curr:
                for (int j = num_cols -1; j >= 0; j--){
                    // T1
                    if ((T1_prev[j+1].value >= T2_prev[j+1].value) && (T1_prev[j+1].value >= T3_prev[j+1].value)){
                        T1_curr[j].value = f(A[i], B[j]) + T1_prev[j+1].value;
                        T1_curr[j].r_row = T1_prev[j+1].r_row;
                        T1_curr[j].r_type = T1_prev[j+1].r_type;

                    } else  if (T2_prev[j+1].value >= T3_prev[j+1].value){
                        T1_curr[j].value = f(A[i], B[j]) + T2_prev[j+1].value;
                        T1_curr[j].r_row = T2_prev[j+1].r_row;
                        T1_curr[j].r_type = T2_prev[j+1].r_type;

                    } else {
                        T1_curr[j].value = f(A[i], B[j]) + T3_prev[j+1].value;
                        T1_curr[j].r_row = T3_prev[j+1].r_row;
                        T1_curr[j].r_type = T3_prev[j+1].r_type;   
                    }
                    // T3
                    int t1 = T1_prev[j].value - (g+h);
                    int t2 = T2_prev[j].value - (g+h);
                    int t3 = T3_prev[j].value - g;
                    if ((t1 >= t2) && (t1 >= t3)){
                        T3_curr[j].value = t1;
                        T3_curr[j].r_row = T1_prev[j].r_row;
                        T3_curr[j].r_type = T1_prev[j].r_type;

                    } else  if (t2 >= t3){
                        T3_curr[j].value = t2;
                        T3_curr[j].r_row = T2_prev[j].r_row;
                        T3_curr[j].r_type = T2_prev[j].r_type;

                    } else {
                        T3_curr[j].value = t3;
                        T3_curr[j].r_row = T3_prev[j].r_row;
                        T3_curr[j].r_type = T3_prev[j].r_type;   
                    }

                    // T2
                    t1 = T1_curr[j+1].value - (g+h);
                    t2 = T2_curr[j+1].value - g;
                    t3 = T3_curr[j+1].value - (g+h);
                    if ((t1 >= t2) && (t1 >= t3)){
                        T2_curr[j].value = t1;
                        T2_curr[j].r_row = T1_curr[j+1].r_row;
                        T2_curr[j].r_type = T1_curr[j+1].r_type;

                    } else if (t2 >= t3){
                        T2_curr[j].value = t2;
                        T2_curr[j].r_row = T2_curr[j+1].r_row;
                        T2_curr[j].r_type = T2_curr[j+1].r_type;

                    } else {
                        T2_curr[j].value = t3;
                        T2_curr[j].r_row = T3_curr[j+1].r_row;
                        T2_curr[j].r_type = T3_curr[j+1].r_type;   
                    }
                    if (j == num_cols -1){
                        T1_curr[num_cols -1].r_row = i;
                        T1_curr[num_cols -1].r_type = 1;
                        T2_curr[num_cols -1].r_row = i;
                        T2_curr[num_cols -1].r_type = 2;
                        T3_curr[num_cols -1].r_row = i;
                        T3_curr[num_cols -1].r_type = 3;
                    }
                }

                // Share last value:
                if (id != col_k1+1){ // Not last row
                    std::unique_lock<std::mutex> lock(working_mutexes[id+2]);
                    while (working[id+2] != 1){
                        update.wait(lock);
                    }
                    (*sharingT)[0][id-2] = T1_curr[0];
                    (*sharingT)[1][id-2] = T2_curr[0];
                    (*sharingT)[2][id-2] = T3_curr[0];
                    {
                        std::lock_guard<std::mutex> lock(working_mutexes[id-2]);
                        working[id-2] = 2;
                    }
                    update.notify_all();
                }

                std::swap(T1_prev, T1_curr);
                std::swap(T2_prev, T2_curr);
                std::swap(T3_prev, T3_curr);

            }
        }

        {
            std::lock_guard<std::mutex> lock(working_mutexes[id]);
            working[id] = 4;
        }
        update.notify_all();


        // SYNCHRONISATION POINT: all partial opt's are computed
        if (ishead){    // Current thread is Head
            std::unique_lock<std::mutex> lock(working_mutexes[id]);
            while (*std::min_element(working.begin() + col_k1*2, working.begin() + col_k2*2) != 4) {
                update.wait(lock);
            }

            // Get Max OPT
            cell max_opt(INF);
            int leftmost_col;

            for (int j = 0; j < num_blocks; j++){
                if ((*sharingOpt)[j].value > max_opt.value){
                    max_opt = (*sharingOpt)[j];
                    leftmost_col = j;
                }
            }

            // SPLIT PROBLEM
            int r1 = max_opt.origin_row  + row_k1;
            int r2 = max_opt.r_row + row_k1; // rows where we split the problem [get using next and prev]
            int t1 = max_opt.origin_type;
            int t2 = max_opt.r_type; // types where we split the problem

            for (int thread_id = col_k1; thread_id < leftmost_col + 1; thread_id++){
                info[thread_id] = Info(col_k1, leftmost_col, row_k1, r1, s, t1);
            }
            for (int thread_id = leftmost_col +1 ; thread_id < col_k2 + 1; thread_id++){
                info[thread_id] = Info(leftmost_col + 1,col_k2, r2, row_k2, t2, e);
            }
            info[leftmost_col] = Info(leftmost_col, leftmost_col, r1, r2, t1, t2);
            info[leftmost_col+1] = Info(leftmost_col, leftmost_col, r1, r2, t1, t2);

            for (int k = col_k1; k < col_k2+1; k++){
                {
                    std::lock_guard<std::mutex> lock(working_mutexes[2*col_k1]);
                    working[2*col_k1] = 1;
                }
                {
                    std::lock_guard<std::mutex> lock(working_mutexes[2*col_k1+1]);
                    working[2*col_k1+1] = 1;
                }
            }
            update.notify_all();
        } else {
            std::unique_lock<std::mutex> lock(working_mutexes[id+2]);
            while (working[id] == 4){
                update.wait(lock);
            }
        }
        // UPDATE DATA:
        col_k1 = info[id].first; // index of first special column (in group)
        col_k2 = info[id].last; // index of last special column (in group)
        row_k1 = info[id].top_row; // index of first special row (in group)
        row_k2 = info[id].bottom_row; // index of last special row (in group)
        s = info[id].s;
        e = info[id].e;
    }
    if (id%2 == 1){ // Even threads solve subproblem
        return;
    }

    // SOLVE SUBPROBLEM


}

int run(std::vector<char> A, std::vector<char> B, int p){
    /*
    We want p even
    */
    if (p % 2 == 1){
        std::cout << "ERROR: Only accept even number of processors." << std::endl;
        return 1;
    }

    std::vector<std::thread> workers(p);
    std::vector<result> results(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1); // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<int>>* sharingT = new std::vector<std::vector<int>>(3, std::vector<int>(p));
    std::vector<std::vector<int>>* sharingTRev = new std::vector<std::vector<int>>(3, std::vector<int>(B.size()));
    std::vector<cell>* sharingOpt = new std::vector<cell>(p/2, cell(INF));
    std::vector<Info>* info = new std::vector<Info>(p/2, Info(0,p/2, 0, B.size()/p, -1, -1));

    for (int i = 0; i<p; i++){
        workers[i] = std::thread(thread_f, &A, &B, p,i, A.size(), B.size(), std::ref(Working), std::ref(working_mutexes), std::ref(update), sharingTRev, sharingOpt, std::ref(info), std::ref(results[i]));
    }
    int TotalScore = 0;
    alignment FinalAlignment;
    for (int i = 0; i < p; i++){
        workers[i].join();
        TotalScore += results[i].score;
        FinalAlignment += results[i].al;
    }

    // join solutions
}

int main(){
    std::vector<char> A = readFastaSequence("sequences/insulin_bovin.fasta");
    std::vector<char> B = readFastaSequence("sequences/insulin_homo.fasta");
    run(A, B, 2);
    return 0;
}