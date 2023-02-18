#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <thread>
#include <climits>
#include <queue>
#include <deque>
#include <chrono>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>
#include "hungarian.hpp"
#include "history.hpp"
#include "active_tree.hpp"

using namespace std;

extern bool enable_threadstop;

struct int_64 {
    int val = 0;
    bool explored = false;
    bool padding[59];
};

struct bool_64 {
    bool val = false;
    bool padding[63];
};

struct unsigned_long_64 {
    unsigned long val = 0;
    bool padding[56];
};

struct mutex_64 {
    mutex lck;
    bool padding[24];
};

enum restart_criteria {
    BEST_SOLUTION_BASED,
    HISTORY_UPDATE_BASED,
};

enum pool_dest {
    GLOBAL_POOL,
    LOCAL_POOL,
};

enum space_state {
    NOT_PROMISING,
    IS_PROMISING,
};

enum space_ranking {
    ABANDONED,
    PROMISING,
    UNEXPLORED,
};

enum abandon_action {
    EXPLOITATION_PROMISE,
    EXPLORE_UNKNOWN,
    CONTINUE_EXPLORATION,
};

enum transformation_option {
    T_ENABLE,
    T_DISABLE,
};

struct request_packet {
    int target_lastnode;
    int target_depth;
    int target_prefix_cost;
    int target_thread;
    boost::dynamic_bitset<> key;
};

class load_stats {
    public:
        bool out_of_work;
        bool stop_checked;
        bool resume_checked;
        bool stop_sig;
        vector<int>* cur_sol;
        space_state space_ranking;
        abandon_action next_action;
        std::chrono::time_point<std::chrono::system_clock> found_time;
        int data_cnt;

        bool padding[44];
        
        load_stats() {
            cur_sol = NULL;
            stop_checked = false;
            resume_checked = false;
            out_of_work = false;
            stop_sig = false;
            next_action = abandon_action::EXPLORE_UNKNOWN;
            data_cnt = 0;
        }
};

class edge {
    public:
        int src;
        int dest;
        int weight;
        edge(int x, int y, int z): src{x},dest{y},weight{z} {}
        bool operator==(const edge &rhs) const {return this->src == rhs.src;}
};

class node {
    public:
        int n = -1;
        int lb = -1;
        int nc = -1;
        int partial_cost = -1;
        bool node_invalid = false;
        HistoryNode* his_entry = NULL;
        Active_Node* act_entry = NULL;
        bool pushed = false;;
        node(int x, int y): n{x},lb{y} {}
};

struct promise_stats {
    int lb = INT_MAX;
    unsigned depth = 0;
};

struct speed_resinfo {
    int tid = -1;
    int u = -1;
    int v = -1;
    int start = -1;
    int finish = -1;
    bool compute = false;
    deque<node>* lbproc_list = NULL;
    bool padding[31];
};

template <typename T>
class atomwrapper
{
    private:
        atomic<T> data;
    public:
        atomwrapper() {data.store(false);}
        atomwrapper(const atomic<T> &val) {data.store(val.load());}
        T load() {return data.load();};
        void store(T val) {data.store(val);}
};

class instrct_node {
    public:
        vector<int> sequence;
        int originate = -1;
        int load_info = -1;
        int parent_lv = -1;
        bool* invalid_ptr = NULL;
        bool deprecated = false;

        Active_Path partial_active_path;
        HistoryNode* root_his_node;
        instrct_node(vector<int> sequence_src,int originate_src, int load_info_src, 
                     int best_costrecord_src, Active_Path temp_active_path, 
                     HistoryNode* temp_root_hisnode) {
            sequence = sequence_src;
            originate = originate_src;
            load_info = load_info_src;
            parent_lv = best_costrecord_src;
            partial_active_path = temp_active_path;
            root_his_node = temp_root_hisnode;
        }
};

struct lptr_64 {
    deque<instrct_node> *local_pool = NULL;
    bool padding[56];
};

class sub_pool {
    public:
        atomic<bool>* enable;
        deque<deque<instrct_node>> space;
        sub_pool() {
            enable = new atomic<bool>();
            *enable = false;
        }
        bool is_empty() {
            int size = 0;
            for (auto sub_space : space) {
                size += sub_space.size();
            }

            if (size > 0) return false;
            else return true;

            return false;
        }
};

struct sop_state {
    vector<int> depCnt;
    vector<int> taken_arr;
    vector<int> cur_solution;
    vector<HistoryNode*> history_arr;
    Hungarian hungarian_solver;
    std::chrono::time_point<std::chrono::system_clock> interval_start;
    pair<boost::dynamic_bitset<>,int> key;
    HistoryNode* cur_parent_hisnode = NULL;
    int cur_cost = 0;
    int initial_depth = 0;
    int suffix_cost = 0;
    ////////////////////
    int originate = -1;
    int load_info = -2;
    ///////////////////
};

struct recur_state {
    sop_state state;
    Active_Path atree;
};

struct Global_Pool {
    sub_pool Promise;
    deque<instrct_node> Abandoned;
    deque<instrct_node> Unknown;
};

class solver {
    private:        
        deque<instrct_node> wrksteal_pool;
        deque<instrct_node> *local_pool = NULL;
        vector<recur_state> recur_stack;
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor(vector<int>* partial_solution);
        
        Active_Allocator Allocator;
        Active_Path cur_active_tree;
        bool abandon_work = false;
        bool abandon_share = false;
        bool grabbed = false;
        int node_count = 0;
        int restart_group_id = -1;
        int mg_id = -1;
        bool speed_search = false;
        int lb_curlv = INT_MAX;

        sop_state problem_state;
        sop_state back_up_state;
        HistoryNode* current_hisnode;

        int thread_id = -1;

        //Restart
        int concentrate_lv = 0;

        //Thread Stopping
        int stop_depth = -1;
        int last_node = -1;

        bool stop_init = false;

        bool Wlkload_Request();
        bool HistoryUtilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,HistoryNode** entry, int cost);
        bool check_satisfiablity(int* local_cost, vector<int>* tour);
        bool Split_local_pool();
        bool push_to_pool(pool_dest decision, space_ranking problem_property);
        bool Split_level_check(deque<sop_state>* solver_container);
        bool Is_promisethread(int t_id);
        bool Grab_from_GPQ(bool reserve);
        bool grab_restartnode();
        bool Check_Local_Pool(deque<node>& enumeration_list,deque<node>& curlocal_nodes);
        bool assign_workload(node& transfer_wlkload, pool_dest destination, space_ranking problem_property, HistoryNode* temp_hisnode);
        bool compare_sequence(vector<int>& sequence, int& target_depth);
        bool Steal_Workload();
        bool EnumerationList_PreProcess(deque<node>& enumeration_list,deque<node>& curlocal_nodes);
        int get_maxedgeweight();
        int dynamic_hungarian(int src, int dest);
        int shared_enumerate(int i);
        void enumerate();
        int get_mgid();
        void transfer_wlkload();
        void retrieve_from_local(deque<node>& curlocal_nodes, deque<node>& enumeration_list);
        void initialize_node_sharing();
        void Thread_Selection();
        void Shared_Thread_Selection();
        void Generate_SolverState(instrct_node& sequence_node);
        void push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,HistoryNode** entry,bool backtracked);
        void check_workload_request(int i);
        void notify_finished();
        size_t transitive_closure(vector<vector<int>>& isucc_graph);
        void solve_parallel(int thread_num, int pool_size);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
        void Select_GroupMember(int unpromise_cnt);
        void Check_Restart_Status(deque<node>& enumeration_list, deque<node>& curlocal_nodes);
        void ThreadStopDelete_Pool(int& target_depth);
        void StopCurThread(int target_depth);
        void Local_PoolConfig(float precedence_density);
        void regenerate_hungstate();
        void CheckStop_Request();
    public:
        void solve(string filename,int thread_num);
        void assign_parameter(vector<string> setting);
};

#endif
