#include "alg.h"
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <limits>
#include "read_fasta.cpp"
#include <condition_variable>
#include <cmath>

int INF = std::numeric_limits<int>::min();
int SUP = std::numeric_limits<int>::max();

const int h = 2; // gap_creation_penalty
const int g = 1;  // gap_penalty


class cell{
public:
    cell(int v, int n, int p){value = v; next = n; prev = p;}
    int value;
    int next, prev;
};
    
class result{
    int TO_DO;
};

class Info{
public:
    Info(int f, int l, int t, int b, int st, int et) {
        first = f;
        last = l;
        top_row = t;
        bottom_row = b;
        s = st;
        e = et;
    }
    int first, last, top_row, bottom_row, s, e;
};

int f(char a, char b){
    if(a == b){return 2;} 
    return 0;
};

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
    std::vector<std::mutex>& working_mutexes, std::condition_variable& update, std::vector<std::vector<int>>* sharingT,
    std::vector<std::vector<int>>* sharingTRev, std::vector<cell>* sharingOpt, std::vector<Info>& info,
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
        std::vector<int> T1_curr(num_cols + 1);
        std::vector<int> T2_curr(num_cols + 1);
        std::vector<int> T3_curr(num_cols + 1);

        std::vector<int> T1_prev(num_cols + 1, INF);
        std::vector<int> T2_prev(num_cols + 1, INF);
        std::vector<int> T3_prev(num_cols + 1, INF);

        if(id % 2 == 0){
            int m = num_cols; // amount of columns per thread

            // INITIALIZE T_prev for all threads:
            for (int j = 0; j < num_cols; j++){          
                    if (!((s == 1) || (s == 3))){ T2_prev[j] = -h-g*(j + 1);}
                    if (s == -2){ T2_prev[j] += h;}
            }

            if(ishead){
                // [i'-1, j'-1]
                if (s == -1){ T1_prev[0] = 0; }
                if (s == -2){ T2_prev[0] = 0;}
                if (s == -3){ T3_prev[0] = 0;}
                // [i, j'-1]
                
            }

            for (int i = 0; i < num_rows; i++){ // COMPUTE ALL ROWS
                // Get T_curr[0]:
                if (ishead){
                    T1_curr[0] = INF; // NEEDED?
                    T2_curr[0] = INF; // NEEDED?
                    T3_curr[0] = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[0] = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[0] += h;}
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
                    T1_curr[j] = f(A[i], B[j]) + std::max({T1_prev[j-1], T2_prev[j-1], T3_prev[j-1]});
                    T3_curr[j] = std::max({T1_prev[j] - (g+h), T2_prev[j] - (g+h), T3_prev[j] - g});
                    T2_curr[j] = std::max({T1_curr[j-1] - (g+h), T2_curr[j-1] - g, T3_curr[j-1] - (g+h)});
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
            cell max_opt(INF,0,0);
            int start = p * id * 2;
            int end = start + num_cols;
            std::vector<int> T1_Rev(sharingTRev->at(0).begin() + start, sharingTRev->at(0).begin() + end);
            std::vector<int> T2_Rev(sharingTRev->at(1).begin() + start, sharingTRev->at(1).begin() + end);
            std::vector<int> T3_Rev(sharingTRev->at(2).begin() + start, sharingTRev->at(2).begin() + end);


            for (int j = 0; j < num_cols; j++){
                int tmp;
                int TRevmax = std::max({T1_Rev[j], T2_Rev[j], T3_Rev[j]});
                int Tmax = std::max({T1_curr[j], T2_curr[j], T3_curr[j]});
                tmp = std::max({Tmax + TRevmax, T2_curr[j] + T2_Rev[j], T3_curr[j] + T3_Rev[j]});
                if (tmp > max_opt.value){
                    max_opt.value = tmp;
                    max_opt.next = j;   // TO DO
                    max_opt.prev = j;   // TO DO
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
                    if (!((e == 1) || (e == 3))){ T2_prev[j] = -h-g*(j + 1);}
                    if (e == -2){T2_prev[j] += h;}
            }

            if(id == col_k2+1){ // Is "sub-head"
                // [i'-1, j'-1]
                if (s == -1){ T1_prev[0] = 0; }
                if (s == -2){ T2_prev[0] = 0;}
                if (s == -3){ T3_prev[0] = 0;}
                // [i, j'-1]
                
            }

            for (int i = num_rows -1; i >= 0; i--){ // COMPUTE ALL ROWS
                // Get T_curr[0]:
                if (id == p){
                    T1_curr[0] = INF; // NEEDED?
                    T2_curr[0] = INF; // NEEDED?
                    T3_curr[0] = INF;
                    if (!((s == 1) || (s == 2))){ T3_curr[0] = -h-g*(i + 1);}
                    if (s == -3){ T3_curr[0] += h;}
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
                    T1_curr[j] = f(A[i], B[j]) + std::max({T1_prev[j+1], T2_prev[j+1], T3_prev[j+1]});
                    T3_curr[j] = std::max({T1_prev[j] - (g+h), T2_prev[j] - (g+h), T3_prev[j] - g});
                    T2_curr[j] = std::max({T1_curr[j+1] - (g+h), T2_curr[j+1] - g, T3_curr[j+1] - (g+h)});
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
            cell max_opt(INF,0,0);
            int leftmost_col;

            for (int j = 0; j < num_blocks; j++){
                if ((*sharingOpt)[j].value > max_opt.value){
                    max_opt.value = (*sharingOpt)[j].value;
                    max_opt.next = j;   // TO DO
                    max_opt.prev = j;   // TO DO
                    leftmost_col = j;
                }
            }

            // SPLIT PROBLEM
            int r1, r2; // rows where we split the problem [get using next and prev]
            int t1, t2; // types where we split the problem


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

int main(A, B, p){
    /*
    We want p even
    */
    if (p % 2 == 1){
        std::cout << "ERROR: Only accept even number of processors." << std::endl;
        return 1;
    }

    std::vector<std::thread> workers(num_threads);
    std::vector<result> results(num_threads);
    std::condition_variable update; // notifies if there's a change to Working
    std::vector<int> Working(p, 1); // vector of size p, Working[i] = 1 if thread i is working, 0 otherwise
    std::vector<std::mutex> working_mutexes(p);
    std::vector<std::vector<int>>* sharingT = new std::vector<std::vector<int>>(3, std::vector<int>(p));
    std::vector<std::vector<int>>* sharingTRev = new std::vector<std::vector<int>>(3, std::vector<int>(B.size()));
    std::vector<cell>* sharingOpt = new std::vector<cell>(p/2, cell(INF,0,0));
    std::vector<Info>* info = new std::vector<Info>(p/2, Info(0,p/2, 0, B.size()/p, -1, -1));

    for (int i = 0; i<p; i++){
        workers[i] = std::thread(thread_f, &A, &B, p,i, A.size(), B.size(), std::ref(Working), std::ref(working_mutexes), std::ref(update), sharingTRev, sharingOpt, std::ref(info), std::ref(results[i]));
    }

    for (int i = 0; i < p; i++){
        workers[i].join();
    }

    // join solutions
}
