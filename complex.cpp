#include "complex.h"
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

/* Function does not work correctly */
result solve_subproblem(extended_P sub_problem, int col_offset, int row_offset, int overlap){ // one or two sub_problems

    str A = sub_problem.A;
    str B = sub_problem.B;
    int start_type = sub_problem.s;
    int end_type = sub_problem.e;

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
        if (end_type == 1){ res.score = T1[n-overlap][m];}
        else if (end_type == 2) {res.score = T2[n-overlap][m];}
        else {res.score = T3[n-overlap][m];}
    } else {
        res.score = std::max({ T1[n-overlap][m], T2[n-overlap][m], T3[n-overlap][m] });
    } 

    /* doing the traceback */
    int state = from(res.score, T1[n][m], T2[n][m], T3[n][m]);
    int i(n), j(m);

    i-=overlap;

    while (i > 0 && j > 0) {
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

    while(i>0){
        res.al.data.push_back({i-1 + row_offset, -1});
        i--;
    }
    while(j>0){
        res.al.data.push_back({-1, j-1 + col_offset});
        j--;
    }

    std::reverse(res.al.data.begin(), res.al.data.end());
    return res;
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
    sharingTRev : to share values of TRev to even processors
    sharingOpt : to share opt values between processors
    Info : what subproblem each thread is working on
*/
void thread_f(const char* A_ptr, const char* B_ptr, int p, int id, int n, int m, std::vector<int>& working,
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<cell>>* sharingT,
    std::vector<std::vector<cell>> &sharingTRev, std::vector<cell>* sharingOpt, std::vector<Info>& info,
    result& result){
    int leftmost_col;
    int iter = 0;
    int num_blocks = p / 2;
    if (p==1){num_blocks=1;}

    int B_len;
    if ((id == p-2) || (id == p-1)){
        B_len = (m / num_blocks) + (m % num_blocks);

    } else {
        B_len = m/num_blocks;
    }

    int col_k1 = info[id].first; // index of first special column (in group)
    int col_k2 = info[id].last; // index of last special column (in group)
    int row_k1 = info[id].top_row; // index of first row (in group)
    int row_k2 = info[id].bottom_row; // index of last row (in group)
    int s = info[id].s;
    int e = info[id].e;

    while (((col_k2 - col_k1 > 1) && (row_k1 +1 != row_k2)) && ((p > 2)&& (B_len >10))){
        // Decomposing

        int row_midk = ceil((row_k1 + row_k2)/2);
        bool ishead = (2*col_k1 == id);
        B_len += info[id].plus;
        std::vector<char> B(B_ptr + (m/num_blocks)*(id/2), B_ptr + (m/num_blocks)*(id/2) + B_len); // Copy B part

        // Initialize curr and prev
        std::vector<cell> T1_curr(B_len + 1, INF);
        std::vector<cell> T2_curr(B_len + 1, INF);
        std::vector<cell> T3_curr(B_len + 1, INF);

        std::vector<cell> T1_prev(B_len + 1, INF);
        std::vector<cell> T2_prev(B_len + 1, INF);
        std::vector<cell> T3_prev(B_len + 1, INF);

        if(id % 2 == 0){

            // Copy Section of A we want
            int num_rows = (row_midk - row_k1);
            std::vector<char> A(A_ptr + row_k1, A_ptr + row_k1 + num_rows); // Copy A part

            // INITIALIZE T_prev for all threads:
            for (int j = 0; j < B_len; j++){          
                    if (!((s == 1) || (s == 3))){ T2_prev[j].value = -h-g*(j + 1);}
                    if (s == -2){ T2_prev[j].value += h;}
            }

            if(ishead){
                if (s == -1){ T1_prev[0].value = 0; }
                if (s == -2){ T2_prev[0].value = 0;}
                if (s == -3){ T3_prev[0].value = 0;}
            }

            for (int i = 0; i < num_rows; i++){ // COMPUTE ALL ROWS
                // Get T_curr[0]:
                if (ishead){
                    T1_curr[0].value = INF;
                    T2_curr[0].value = INF;
                    T3_curr[0].value = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[0].value = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[0].value += h;}
                } else {
                    { // Read from shared
                        std::unique_lock<std::mutex> lock(working_mutexes[id]);
                        while (working[id] != 2){
                            update.wait(lock);
                        }
                        T1_curr[0] = (*sharingT)[0][id];
                        T2_curr[0] = (*sharingT)[1][id];
                        T3_curr[0] = (*sharingT)[2][id];
                        working[id] = 1;
                    }
                    update.notify_all();
                }

                // Compute T_curr:
                for(int j = 1; j<B_len; j++){
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
                if (id != col_k2*2){ // If not last column in group
                    {
                        std::unique_lock<std::mutex> lock(working_mutexes[id+2]);
                        while (working[id+2] != 1){
                            update.wait(lock);
                        }
                        (*sharingT)[0][id+2] = T1_curr[B_len-1];
                        (*sharingT)[1][id+2] = T2_curr[B_len-1];
                        (*sharingT)[2][id+2] = T3_curr[B_len-1];
                        working[id+2] = 2;
                    }
                    update.notify_all();
                }

                T1_prev = T1_curr;
                T2_prev = T2_curr;
                T3_prev = T3_curr;
            }
            // T is fully computed
            {
                // Wait until TRev Computed for this set of columns
                std::unique_lock<std::mutex> lock(working_mutexes[id+1]);
                while (working[id+1] != 4){
                    update.wait(lock);
                }
            }
            // Compute opt
            cell max_opt(INF);

            std::vector<cell> T1_Rev(B_len);
            std::vector<cell> T2_Rev(B_len);
            std::vector<cell> T3_Rev(B_len);

            { // Get TRev values for this block
                std::unique_lock<std::mutex> lk(working_mutexes[0]);
                auto& row0 = sharingTRev[0];
                auto& row1 = sharingTRev[1];
                auto& row2 = sharingTRev[2];

                int start_idx = (m / num_blocks) * (id / 2);

                for(int j = 0; j<B_len; j++){
                   T1_Rev[j] = sharingTRev[0][start_idx + j];
                   T2_Rev[j] = sharingTRev[1][start_idx + j];
                   T3_Rev[j] = sharingTRev[2][start_idx + j];
                }
            }

            // Compute OPT
            for (int j = 0; j < B_len; j++){
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
                    o_row = T1_curr[j].origin_row;
                    o_type = T1_curr[j].origin_type;
                    Tmax = T1_curr[j].value;
                } else if (T2_curr[j].value >= T3_curr[j].value){
                    o_row = T2_curr[j].origin_row;
                    o_type = T2_curr[j].origin_type;
                    Tmax = T2_curr[j].value;
                } else {
                    o_row = T3_curr[j].origin_row;
                    o_type = T3_curr[j].origin_type;
                    Tmax = T3_curr[j].value;
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
            
            {// Share Local Opt
                std::lock_guard<std::mutex> lock(working_mutexes[0]);
                (*sharingOpt)[id/2] = max_opt;
            }
            { // Update Working status
                std::lock_guard<std::mutex> lock(working_mutexes[id]);
                working[id] = 4;
            }
            update.notify_all();


        } else {
            // Copy Section of A we want
            int num_rows = (row_k2 - row_midk);
            std::vector<char> A(A_ptr + row_midk, A_ptr + row_midk + num_rows); // Copy A part

            // INITIALIZE T_prev for all threads:
            for (int j = B_len -1; j >= 0; j--){          
                    if (!((e == 1) || (e == 3))){ T2_prev[j].value = -h-g*(j + 1);}
                    if (e == -2){T2_prev[j].value += h;}
            }

            if(id == (2*col_k2)+1){ // Is "sub-head"
                if (s == -1){ T1_prev[0].value = 0; }
                if (s == -2){ T2_prev[0].value = 0;}
                if (s == -3){ T3_prev[0].value = 0;}
                
            }

            for (int i = num_rows -1; i >= 0; i--){ // COMPUTE ALL ROWS
                // Get TRev_curr[0]:
                if (id == (col_k2*2)+1){ // Rev Head
                    T1_curr[B_len -1].value = INF;
                    T2_curr[B_len -1].value = INF;
                    T3_curr[B_len -1].value = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[B_len -1].value = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[B_len -1].value += h;}
                } else {
                    {
                        std::unique_lock<std::mutex> lock(working_mutexes[id]);
                        while (working[id] != 2){
                            update.wait(lock);
                        }
                        T1_curr[B_len] = (*sharingT)[0][id];
                        T2_curr[B_len] = (*sharingT)[1][id];
                        T3_curr[B_len] = (*sharingT)[2][id];
                        working[id] = 1;
                    }
                    update.notify_all();
                }

                // Compute T_curr:
                for (int j = B_len -1; j >= 0; j--){
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
                    if (j == B_len -1){
                        T1_curr[j].r_row = i;
                        T1_curr[j].r_type = 1;
                        T2_curr[j].r_row = i;
                        T2_curr[j].r_type = 2;
                        T3_curr[j].r_row = i;
                        T3_curr[j].r_type = 3;
                    }
                }

                // Share last value:
                if (id != col_k1*2+1){ // If not last [first] column
                    {
                        std::unique_lock<std::mutex> lock(working_mutexes[id-2]);
                        while (working[id-2] != 1){
                            update.wait(lock);
                        }
                        (*sharingT)[0][id-2] = T1_curr[0];
                        (*sharingT)[1][id-2] = T2_curr[0];
                        (*sharingT)[2][id-2] = T3_curr[0];

                        working[id-2] = 2;
                    }
                    update.notify_all();
                }
                T1_prev = T1_curr;
                T2_prev = T2_curr;
                T3_prev = T3_curr;
            }
            // Finished computing TRev

            { // Share TRev values
                std::unique_lock<std::mutex> lk(working_mutexes[0]);
                int start_idx = (m / num_blocks) * (id / 2);
                for(int j = 0; j<B_len; j++){
                    sharingTRev[0][start_idx + j] = T1_curr[j];
                    sharingTRev[1][start_idx + j] = T2_curr[j];
                    sharingTRev[2][start_idx + j] = T3_curr[j];
                }

            }
            { // Update Working status
                std::lock_guard<std::mutex> lock(working_mutexes[id]);
                working[id] = 4;
            }
            update.notify_all();
        }


        // SYNCHRONISATION POINT: all partial opt's are computed and in SharedOpt
        if (ishead){// Current thread is Head
            {   // Wait until all threads in group finished computing SharedOpt
                std::unique_lock<std::mutex> lock(working_mutexes[id]);
                while (*std::min_element(working.begin() + col_k1*2, working.begin() + col_k2*2+1) != 4) {
                    update.wait(lock);
                }
            }
            // Get Max OPT
            cell max_opt(INF);
            for (int j = col_k1; j < col_k2+1; j++){
                // Access shared OPT
                std::lock_guard<std::mutex> lock(working_mutexes[0]);
                if ((*sharingOpt)[j].value > max_opt.value){
                    max_opt = (*sharingOpt)[j];
                    leftmost_col = j;   
                }
            }

            // SPLIT PROBLEM

            int r1 = max_opt.origin_row  + row_k1-1;
            int r2 = max_opt.r_row + row_midk+1; // rows where we split the problem [get using next and prev]
            int t1 = max_opt.origin_type;
            int t2 = max_opt.r_type; // types where we split the problem

            // Update Info for all threads in group:
            for (int thread_id = 2*col_k1; thread_id < 2*leftmost_col; thread_id++){
                info[thread_id] = Info(col_k1, leftmost_col-1, row_k1, r1, s, t1, 1);
            }
            for (int thread_id = 2*leftmost_col +2 ; thread_id < 2*col_k2 + 2; thread_id++){
                info[thread_id] = Info(leftmost_col +1,col_k2, r2, row_k2, t2, e,0);

            }
                info[2*leftmost_col] = Info(leftmost_col, leftmost_col, r1, r2, t1, t2, 1);
                info[2*leftmost_col+1] = Info(leftmost_col, leftmost_col, r1, r2, t1, t2, 1);

            if(leftmost_col == col_k1){
                info[leftmost_col*2] = Info(leftmost_col, leftmost_col, row_k1, r2, t1, t2, 1);
                info[leftmost_col*2+1] = Info(leftmost_col, leftmost_col, row_k1, r2, t1, t2, 1);
            } else if(leftmost_col == col_k2){
                info[leftmost_col*2] = Info(leftmost_col, leftmost_col, r1, row_k2, t1, t2, 1);
                info[leftmost_col*2+1] = Info(leftmost_col, leftmost_col, r1, row_k2, t1, t2, 1);
            }
            // Notify threads to continue working:
            for (int k = col_k1; k < col_k2+1; k++){
                {
                    std::lock_guard<std::mutex> lock(working_mutexes[2*k+1]);
                    working[2*k+1] = 1;
                }
                {
                    std::lock_guard<std::mutex> lock(working_mutexes[2*k]);
                    working[2*k] = 1;
                }
            }
            update.notify_all();
        } else {
            // If not head, wait until info is updated
            std::unique_lock<std::mutex> lock(working_mutexes[id]);
            while (working[id] == 4){
                update.wait(lock);
            }
        }

        // UPDATE PROBLEM DATA:
        col_k1 = info[id].first; // index of first special column (in group)
        col_k2 = info[id].last; // index of last special column (in group)
        row_k1 = info[id].top_row; // index of first special row (in group)
        row_k2 = info[id].bottom_row; // index of last special row (in group)
        s = info[id].s;
        e = info[id].e;
    }

    // SOLVING SUBPROBLEM

    if ((id%2 == 1)||((!(2*col_k1 == id)) || (row_k2-row_k1 < 0))){ // If not head of group - Terminate
        return;
    }

    // Get A and B substrings:
    B_len = info[id].plus;
    if (col_k1 +1 == col_k2){B_len += m/num_blocks;}
    if (2*col_k2 == p-2){
        B_len += (m / num_blocks) + (m % num_blocks);
    } else {
        B_len += m/num_blocks;
    }

    int start_idx = (m / num_blocks) * col_k1;
    
    std::vector<char> B2(B_ptr+start_idx, B_ptr +start_idx+ B_len); // Copy B part

    int num_rows = (row_k2 - row_k1) + info[id].plus;
    std::vector<char> A(A_ptr + row_k1, A_ptr + row_k1 + num_rows); // Copy A part

    extended_P sub_problem;
    sub_problem.A = A;
    sub_problem.B = B2;
    sub_problem.s = s;
    sub_problem.e = e;

    // Solve subproblem:
    result = solve_subproblem(sub_problem, start_idx, row_k1, info[id].plus);
}

int run(std::vector<char> A, std::vector<char> B, int p){

    if (p%2 != 0){p = std::max({1, p-1});} // Algorithm works for even number of threads (since we compute Trev Separately)

    while(B.size()/p <= 10){ // We don't want to have blocks of columns which are too small
        p -= 2;
    }

    // Initialise Variables / Vectors
    std::vector<std::thread> workers(p);
    std::vector<result> results(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1);
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<cell>>* sharingT = new std::vector<std::vector<cell>>(3, std::vector<cell>(p));
    std::vector<std::vector<cell>> sharingTRev(3, std::vector<cell>(B.size()));
    std::vector<cell>* sharingOpt = new std::vector<cell>(p/2, cell(INF));
    std::vector<Info>* info = new std::vector<Info>(p, Info(0,p/2 -1, 0, A.size(), -1, -1,0));

    // Create Threads
    for (int i = 0; i<p-1; i++){
        workers[i] = std::thread(thread_f, A.data(), B.data(), p,i, A.size(), B.size(), std::ref(Working),
        std::ref(working_mutexes), std::ref(update), sharingT, std::ref(sharingTRev), sharingOpt, std::ref(*info), std::ref(results[i]));
    }
    thread_f(A.data(), B.data(), p,p-1, A.size(), B.size(), std::ref(Working),
    std::ref(working_mutexes), std::ref(update), sharingT, std::ref(sharingTRev), sharingOpt, std::ref(*info), std::ref(results[p-1]));

    // Join Workers
    for (int i = 0; i < p-1; i++){
        workers[i].join();
    }

    // Join Partial Solutions:
    int TotalScore = 0;
    alignment FinalAlignment;
    for (int i = 0; i < p; i++){
        TotalScore += results[i].score;
        FinalAlignment += results[i].al;
    }

    std::cout << "Score: " << TotalScore << std::endl;
    output_alignement(FinalAlignment, A, B);
    delete sharingT;
    delete sharingOpt;
    delete info;
    return 0;
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

    std::vector<char> A = readFastaSequence(folder + sequence1_filename);
    std::vector<char> B = readFastaSequence(folder+sequence2_filename);

    // std::cout<< "INPUT SEQUENCES:"<<std::endl;
    // std::cout << "A: " << std::string(A.begin(), A.end()) << "\n"
    //           << "B: " << std::string(B.begin(), B.end()) << std::endl;
    

    // Uncomment to test:

    // for(int k = 1; k < p; k++){
    //     auto start = std::chrono::high_resolution_clock::now();
    //     run(A, B, k);
    //     auto end = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    //     std::cout << "Execution time for " << k << " processors: " << duration.count() << " nanoseconds\n";
    // }

    run(A, B, p);
    return 0;
}