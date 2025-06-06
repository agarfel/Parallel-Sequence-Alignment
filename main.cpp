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

/* function should be correct, watch out for out of bound indices in the trace back maybe but otherwise it works */
result solve_subproblem(extended_P sub_problem, int col_offset, int row_offset, int overlap){ // one or two sub_problems

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
    //if(col_offset!=0){std::cout << 1 << std::endl;}
    // Initialize all cells to -1
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            T1[i][j] = -1;
            T2[i][j] = -1;
            T3[i][j] = -1;
        }
    }

    //if(col_offset!=0){std::cout << 2<< std::endl;}

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
    //if(col_offset!=0){std::cout << 3 << std::endl;}

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
    //if(col_offset!=0){std::cout << 4 << std::endl;}

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
    //if(col_offset!=0){std::cout << 5 << std::endl;}

    if (end_type > 0) {
        if (end_type == 1){ res.score = T1[n-overlap][m];}
        else if (end_type == 2) {res.score = T2[n-overlap][m];}
        else {res.score = T3[n-overlap][m];}
    } else {
        res.score = std::max({ T1[n-overlap][m], T2[n-overlap][m], T3[n-overlap][m] });
    } 
    //if(col_offset!=0){std::cout << 6 << std::endl;}


    /* doing the traceback */
    int state = from(res.score, T1[n][m], T2[n][m], T3[n][m]);
    int i(n), j(m);
    //if(col_offset!=0){std::cout << 7 << std::endl;}
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

        //if(col_offset!=0){std::cout << 8 << std::endl;}


    std::reverse(res.al.data.begin(), res.al.data.end());
    // //if(col_offset!=0){std::cout << 9 << std::endl;}
    // std::cout << "ALIGNMENT LENGTH: " << res.al.data.size() << std::endl;
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

    // Decomposing
    while (((col_k2 - col_k1 > 1) && (row_k1 +1 != row_k2)) && p > 2){
        // {
        // std::lock_guard<std::mutex> lock(working_mutexes[0]);
        // std::cout << "THREAD: " <<id<< " DECOMPOSING" << std::endl;
        // std::cout << col_k2 - col_k1 << std::endl;
        // }
        int row_midk = ceil((row_k1 + row_k2)/2);   //index of middle row
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
            // int m = B_len; // amount of columns per thread
            // if (id == 2*col_k2){
            //     B_len += info[id].plus;
            // }
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
                    T1_curr[0].value = INF; // NEEDED?
                    T2_curr[0].value = INF; // NEEDED?
                    T3_curr[0].value = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[0].value = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[0].value += h;}
                } else {
                    {
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
                // if (id == 4) {

                //         std::lock_guard<std::mutex> lock(working_mutexes[0]);
                //         if (i == 0){
                //                 std::cout << "A: " << std::string(A.begin(), A.end()) << "\n"
                //                 << "B: " << std::string(B.begin(), B.end()) << std::endl;
                //         }
                //         std::cout << "Thread " << id << "\n" << std::endl;
                //         for(int j = 1; j<m; j++){
                //             std::cout << T1_curr[j].value << " ";
                //         }
                // }
                // Share last value:
                if (id != col_k2*2){ // Not last row
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
            // {
            //         std::unique_lock<std::mutex> lock(working_mutexes[0]);

            //         std::cout << "CRY - T of Thread " << id << std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T1_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T2_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T3_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;

            // }
            {
            // Wait until TRev Computed for this set of columns
            std::unique_lock<std::mutex> lock(working_mutexes[id+1]);
            while (working[id+1] != 4){
                update.wait(lock);
            }
            }
            // Copmute opt
            cell max_opt(INF);

            // Safer and clearer: use iterator arithmetic
                std::vector<cell> T1_Rev(B_len);
                std::vector<cell> T2_Rev(B_len);
                std::vector<cell> T3_Rev(B_len);
{
                std::unique_lock<std::mutex> lk(working_mutexes[0]);

                auto& row0 = sharingTRev[0];
                auto& row1 = sharingTRev[1];
                auto& row2 = sharingTRev[2];

                // std::cout << "\n start_idx " << start_idx << " Last: " << start_idx + B_len << " Row0.size() = " << row0.size() << std::endl;
                // std::cout << "Copying thread " << id << ", start_idx = " << start_idx << ", B_len = " << B_len << "\n";
                // std::cout << "T1_curr.size() = " << T1_curr.size() << ", sharingTRev[0].size() = " << sharingTRev[0].size() << std::endl;

                int start_idx = (m / num_blocks) * (id / 2);

                for(int j = 0; j<B_len; j++){
                   T1_Rev[j] = sharingTRev[0][start_idx + j];
                   T2_Rev[j] = sharingTRev[1][start_idx + j];
                   T3_Rev[j] = sharingTRev[2][start_idx + j];

                }
}

            std::vector<int> tmps(B_len);
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
                tmps[j] = tmp;
                if (tmp > max_opt.value){
                    max_opt.value = tmp;
                    max_opt.origin_type = o_type;
                    max_opt.r_type = r_type;
                    max_opt.origin_row = o_row;
                    max_opt.r_row = r_row;
                }
            }
            // Share opt
            {
                std::lock_guard<std::mutex> lock(working_mutexes[0]);
                (*sharingOpt)[id/2] = max_opt;
            }
            {
                std::lock_guard<std::mutex> lock(working_mutexes[id]);
                working[id] = 4;
            }
            update.notify_all();


        } else {
            // Copy Section of A we want
            int num_rows = (row_k2 - row_midk);
            std::vector<char> A(A_ptr + row_midk, A_ptr + row_midk + num_rows); // Copy A part

            // Compute TRev
            // if (id == 2*col_k2 +1){
            //     B_len += info[id].plus;
            // }

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

                    T1_curr[B_len -1].value = INF; // NEEDED?
                    T2_curr[B_len -1].value = INF; // NEEDED?
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
                if (id != col_k1*2+1){ // Not last row
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
            // {
            //         std::unique_lock<std::mutex> lock(working_mutexes[0]);

            //         std::cout << "CRY - TREv of Thread " << id << std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T1_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T2_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;
            //         for(int j = 0; j<B_len; j++){
            //             std::cout << T3_curr[j].value<< " ";
            //         }
            //         std::cout<<std::endl;

            // }
        {
            std::unique_lock<std::mutex> lk(working_mutexes[0]);
            int start_idx = (m / num_blocks) * (id / 2);
            for(int j = 0; j<B_len; j++){
                   sharingTRev[0][start_idx + j] = T1_curr[j];
                   sharingTRev[1][start_idx + j] = T2_curr[j];
                   sharingTRev[2][start_idx + j] = T3_curr[j];
                }
            // std::cout << "CRY" << std::endl;
            // for(int j = 0; j<B_len; j++){
            //     std::cout << sharingTRev[0][start_idx + j].r_type<<
            //        sharingTRev[1][start_idx + j].r_row<<
            //        sharingTRev[2][start_idx + j].r_row << std::endl;
            // }
        }
        {
            std::lock_guard<std::mutex> lock(working_mutexes[id]);
            working[id] = 4;
        }
        update.notify_all();
        }


        // SYNCHRONISATION POINT: all partial opt's are computed
        if (ishead){    // Current thread is Head
            {
                std::unique_lock<std::mutex> lock(working_mutexes[id]);
                while (*std::min_element(working.begin() + col_k1*2, working.begin() + col_k2*2+1) != 4) {
                    update.wait(lock);
                }
            }
            // Get Max OPT
            cell max_opt(INF);

            for (int j = col_k1; j < col_k2+1; j++){

                std::lock_guard<std::mutex> lock(working_mutexes[0]);
                if ((*sharingOpt)[j].value > max_opt.value){
                    max_opt = (*sharingOpt)[j];
                    leftmost_col = j;   
                }
            }
 
                         
            //     {
            //     std::lock_guard<std::mutex> lk(working_mutexes[0]);
            //     std::cout<< "OPT For Columns: " << j << "\n"
            //     << "origin_row: " << (*sharingOpt)[j].origin_row <<"\n"
            //     << "origin_type: " << (*sharingOpt)[j].origin_type <<"\n"
            //     << "r_row: " << (*sharingOpt)[j].r_row <<"\n"
            //     << "r_type: " << (*sharingOpt)[j].r_type <<"\n"
            //     << "value: " << (*sharingOpt)[j].value <<"\n"
            //     << std::endl;
            // }
            
            // {
            //     std::lock_guard<std::mutex> lk(working_mutexes[0]);
            //     std::cout<< "OPT VALUES: " << std::endl;
            //     for (int j= 0; j<num_blocks; j++){
            //         std::cout<<(*sharingOpt)[j].value << " ";
            //     }
            // }
            // SPLIT PROBLEM

            int r1 = max_opt.origin_row  + row_k1-1;
            int r2 = max_opt.r_row + row_midk+1; // rows where we split the problem [get using next and prev]
            int t1 = max_opt.origin_type;
            int t2 = max_opt.r_type; // types where we split the problem
            // {
            //     std::lock_guard<std::mutex> lk(working_mutexes[0]);
            //     std::cout<< "Splitting Problem: " << id << "\n"
            //     << "r1: " << r1 <<"\n"
            //     << "row_k1: " << row_k1 <<"\n"
            //     << "max_opt.origin_row: " << max_opt.origin_row <<"\n"
            //     << "r2: " << r2 <<"\n"
            //     << "row_midk: " << row_midk <<"\n"
            //     << "max_opt.r_row: " << max_opt.r_row <<"\n"
            //     << "t1: " << t1 <<"\n"
            //     << "t2: " << t1 <<"\n"
            //     << "row_k1: " << row_k1 <<"\n"
            //     << "row_k2: " << row_k2 <<"\n"
            //     << "Leftmost: " << leftmost_col <<"\n \n"
            //     << std::endl;
            // }
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
            std::unique_lock<std::mutex> lock(working_mutexes[id]);
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
            // if(id<6){
            //     std::lock_guard<std::mutex> lk(working_mutexes[0]);
            //     std::cout<< "UPDATED INFO - THREAD: " << id << "\n"
            //     << "col_k1: " << col_k1 <<"\n"
            //     << "col_k2: " << col_k2 <<"\n"
            //     << "row_k1: " << row_k1 <<"\n"
            //     << "row_k2: " << row_k2 <<"\n"
            //     << "s: " << s <<"\n"
            //     << "e: " << e <<"\n"
            //     << std::endl;
            // }
            // if(iter == 2){return;}
            // iter++;
    }


    if ((id%2 == 1)||((!(2*col_k1 == id)) || (row_k2-row_k1 < 0))){
        // {
        //     std::lock_guard<std::mutex> lock(working_mutexes[0]);
        //     std::cout << "KILLING THREAD: " <<id<< std::endl;
        // }
        return;
    }
// {
//     std::lock_guard<std::mutex> lock(working_mutexes[0]);
//     std::cout << "THREAD: " <<id<< " SOLVING" << std::endl;
// }
            // {
            //     std::lock_guard<std::mutex> lk(working_mutexes[0]);
            //     std::cout<< "UPDATED INFO - THREAD: " << id << "\n"
            //     << "col_k1: " << col_k1 <<"\n"
            //     << "col_k2: " << col_k2 <<"\n"
            //     << "row_k1: " << row_k1 <<"\n"
            //     << "row_k2: " << row_k2 <<"\n"
            //     << "overlap: " << info[id].plus <<"\n"
            //     << "s: " << s <<"\n"
            //     << "e: " << e <<"\n"
            //     << std::endl;
            // }
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

                {
                std::lock_guard<std::mutex> lk(working_mutexes[0]);
                std::cout<< "SUBPROBLEM INFO - THREAD: " << id << "\n"
                << "A: " << std::string(A.begin(), A.end()) <<"\n"
                << "B: " << std::string(B2.begin(), B2.end()) <<"\n"
                << "s: " << s <<"\n"
                << "e: " << e <<"\n"
                << "row_k1: " << row_k1 <<"\n"
                << "row_k2: " << row_k2 <<"\n"
                << "B_len: " << B_len <<"\n"
                << "Start_idx: " << start_idx <<"\n"
                << std::endl;
            }
    result = solve_subproblem(sub_problem, start_idx, row_k1, info[id].plus);
    // std::cout << "DONE: " << id << std::endl;
            // {
            //  std::lock_guard<std::mutex> lk(working_mutexes[0]);
            // std::vector<char> A_F(A_ptr, A_ptr + n); // Copy B part
            // std::vector<char> B_F(B_ptr, B_ptr + m); // Copy B part
            // output_alignement(result.al, A_F, B_F);

            // }

}

int run(std::vector<char> A, std::vector<char> B, int p){

    if (p%2 != 0){p = std::max({1, p-1});}

    std::vector<std::thread> workers(p);
    std::vector<result> results(p);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1); // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<cell>>* sharingT = new std::vector<std::vector<cell>>(3, std::vector<cell>(p));
    std::vector<std::vector<cell>> sharingTRev(3, std::vector<cell>(B.size()));

    std::vector<cell>* sharingOpt = new std::vector<cell>(p/2, cell(INF));
    std::vector<Info>* info = new std::vector<Info>(p, Info(0,p/2 -1, 0, A.size(), -1, -1,0));

    for (int i = 0; i<p; i++){
        workers[i] = std::thread(thread_f, A.data(), B.data(), p,i, A.size(), B.size(), std::ref(Working),
        std::ref(working_mutexes), std::ref(update), sharingT, std::ref(sharingTRev), sharingOpt, std::ref(*info), std::ref(results[i]));
    }
    int TotalScore = 0;
    alignment FinalAlignment;
    for (int i = 0; i < p; i++){
        workers[i].join();
    }
    for (int i = 0; i < p; i++){
        TotalScore += results[i].score;
        FinalAlignment += results[i].al;
        // output_alignement(results[i].al, A, B);

    }
    std::cout << "Score: " << TotalScore << std::endl;
    output_alignement(FinalAlignment, A, B);
    delete sharingT;
    delete sharingOpt;
    delete info;

    return 0;
}


int main(int argc, char* argv[]) {
    // if (argc != 4) {
    //     std::cerr << "Usage: " << argv[0] << " <sequence1.fasta> <sequence2.fasta> <p>\n";
    //     return 1;
    // }
    std::string folder = "sequences/";

    std::string sequence1_filename = "insulin_homo.fasta"; //argv[1];
    std::string sequence2_filename = "insulin_bovin.fasta";//argv[2];
    int p = 12;//std::stoi(argv[3]);  // Convert string to int
    std::vector<char> A = readFastaSequence(folder + sequence1_filename);
    std::vector<char> B = readFastaSequence(folder+sequence2_filename);
    std::cout<< "INPUT SEQUENCES:"<<std::endl;
    std::cout << "A: " << std::string(A.begin(), A.end()) << "\n"
              << "B: " << std::string(B.begin(), B.end()) << std::endl;
    // for(p = 1; p < 24; p++){
    // auto start = std::chrono::high_resolution_clock::now();
    // run(A, B, p);
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // std::cout << "Execution time for " << p << " processors: " << duration.count() << " nanoseconds\n";
    // }

    // auto start = std::chrono::high_resolution_clock::now();
    run(A, B, 10);
    // run(A, B, 12);

    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // std::cout << "Execution time for " << p << " processors: " << duration.count() << " nanoseconds\n";
    return 0;
}

/*

Score: 169
Length Sequence A : 110
Length Sequence B : 110
Sequence A : MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
Sequence B : MALWTRLRPLLALLALWPPPPARAFVNQHLCGSHLVEALYLVCGERGFFYTPKARREVEGPQVGALELAGGPGAG-----GLEGPPQKRGIVEQCCASVCSLYQLENYCN




*/

