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

/* Basic structure of a specific partial path, a node in the enumeration tree, not including the full information in an instrct_node. */
class node {
    public:
        int n = -1; //node number
        int lb = -1; //lower bound cost of this node
        int nc = -1; //cost to get from the parent to this node
        int partial_cost = -1; //total cost of path to parent
        unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node
        bool node_invalid = false; //investigate
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

/* Entry in local or global pools containing all the information about a path to regenerate an sop_state. */
class instrct_node {
    public:
        vector<int> sequence;
        int originate = -1;
        int load_info = -1;
        int parent_lv = -1;
        bool* invalid_ptr = NULL; //investigate
        bool deprecated = false; //For Thread Stopping. if this node exists in a redundant subspace and so does not need to be processed
        unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node

        Active_Path partial_active_path;
        HistoryNode* root_his_node;
        instrct_node(vector<int> sequence_src,int originate_src, int load_info_src, 
                     int best_costrecord_src, unsigned long long value, Active_Path temp_active_path, 
                     HistoryNode* temp_root_hisnode) {
            sequence = sequence_src;
            originate = originate_src;
            load_info = load_info_src;
            parent_lv = best_costrecord_src;
            current_node_value = value;
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

/* All the information necessary about the current node in the enumeration tree */
struct sop_state {
    vector<int> depCnt;
    vector<int> taken_arr; //a vector of bits for each node in the graph, 1 for taken, 0 for not taken
    vector<int> cur_solution;
    Hungarian hungarian_solver;
    std::chrono::time_point<std::chrono::system_clock> interval_start;
    pair<boost::dynamic_bitset<>,int> key;
    HistoryNode* cur_parent_hisnode = NULL;
    int cur_cost = 0;
    int initial_depth = 0; //depth at which enumeration began once GPQ was initially filled
    int suffix_cost = 0;
    ////////////////////
    int originate = -1;
    int load_info = -2;
    ///////////////////
    unsigned long long current_node_value = -1; //the portion out of ULLONG_MAX of the working tree that is under this node (the partial path represented by this state)
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

/* The local information that each individual thread uses to process paths */
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

        sop_state problem_state; //this thread's current state
        sop_state back_up_state;
        HistoryNode* current_hisnode;

        int thread_id = -1;

        //Restart
        int concentrate_lv = 0;

        //Thread Stopping
        int stop_depth = -1;
        int last_node = -1;

        bool stop_init = false; //whether this thread has ever been stopped before

        /* Build graph based on .sop input file. */
        void retrieve_input(string filename);
        /* Removes redundant edges from the cost graph. */
        void transitive_redundantcy();
        /* Computes the transitive closure of the graph. Used for calculating precedence density. */
        size_t transitive_closure(vector<vector<int>>& isucc_graph);
        /* Initialize local pools. */
        void Local_PoolConfig(float precedence_density);
        /* Sort the cost graph in ascending order. Required for nearest neighbor heuristic. */
        void sort_weight(vector<vector<edge>>& graph);
        /* Find the highest edge weight in the entire cost graph. Required for Hungarian algorithm. */
        int get_maxedgeweight();
        bool check_satisfiablity(int* local_cost, vector<int>* tour);
        bool Split_local_pool();

        /* Called from solver::solve, divides work among global pool and each thread, and begins the threads with calls to solver::enumerate. */
        void solve_parallel(int thread_num, int pool_size); 
        /* Returns true if any solver has a depth different than any other, false otherwise. */
        bool Split_level_check(deque<sop_state>* solver_container);

        /* Recursive function that each thread runs to process its assigned sections of the enumeration tree. */
        void enumerate();
        void notify_finished();
        int shared_enumerate(int i);
        /* Checks whether the next item to be enumerated should instead be pruned. Returns true if the node was pruned, false otherwise. */
        bool EnumerationList_PreProcess(deque<node>& enumeration_list,deque<node>& curlocal_nodes);
        /* Compute the lower bound based on the current hung state with this node added*/
        int dynamic_hungarian(int src, int dest);
        /* Compares the current node to an entry in the history table, if possible. Returns false if this node should be pruned, true otherwise. */
        bool HistoryUtilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,HistoryNode** entry, int cost);
        /* Add an entry to the history table. */
        void push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,HistoryNode** entry,bool backtracked);

        /* Build an sop_state based off the information in the given node. */
        void Generate_SolverState(instrct_node& sequence_node);
        /* Build a hungarian solver state based upon the problem_state. Used in Generate_SolverState. */
        void regenerate_hungstate();

        /* Moves a node into a global or local pool. Returns true if the node was added, false if the node was pruned. 
            Progress tracking for pruned nodes is done internally. */
        bool assign_workload(node& transfer_wlkload, pool_dest destination, space_ranking problem_property, HistoryNode* temp_hisnode);
        /* Push all nodes from this thread's local pool into a global pool. */
        bool push_to_pool(pool_dest decision, space_ranking problem_property);

        //Thread Restart
            void Check_Restart_Status(deque<node>& enumeration_list, deque<node>& curlocal_nodes);
            int get_mgid();
            /* Categorize threads by space_state an issue thread restart. Then calls Select_GroupMember. */
            void Thread_Selection();
            void Shared_Thread_Selection();
            /* Determines if a thread is promising. Returns true if the thread is promising and unexplored, false otherwise. */
            bool Is_promisethread(int t_id);
            /* Decide upon abandon_action for each thread. */
            void Select_GroupMember(int unpromise_cnt);

        //Work Stealing
            /* If the local pool is too small, takes some nodes from enumeration list and moves them into local pool before continuing down the tree. */
            bool Check_Local_Pool(deque<node>& enumeration_list,deque<node>& curlocal_nodes);
            /* Undoes solver::Check_Local_Pool. When backtracking, take back nodes from local pool into the enumeration list in order to ensure proper search order. */
            void retrieve_from_local(deque<node>& curlocal_nodes, deque<node>& enumeration_list);
            /* Check if you should find more work before backtracking. Returns false if work was found and a new state generated to be followed before backtracking. */
            bool Wlkload_Request();
            void initialize_node_sharing();
            bool grab_restartnode();
            /* Try to grab a node from global pool. Returns false if work was found and a new state generated, true otherwise. */
            bool Grab_from_GPQ(bool reserve);
            /* Repeatedly calls transfer_wlkload in order to find work for an idle thread. */
            bool Steal_Workload();
            /* Idle thread randomly chooses a victim thread and takes some nodes from its local_pool to the wrksteal_pool. */
            void transfer_wlkload();
            //void check_workload_request(int i);

        //Thread Stopping
            /* For Thread Stopping. Check the request buffer. */
            void CheckStop_Request();
            /* Process a request and determine if the current thread needs to be stopped. */
            bool thread_stop_check(int target_prefix_cost, int target_depth, int target_lastnode, boost::dynamic_bitset<>& src_key);
            /* Stop the current thread. */
            void StopCurThread(int target_depth);
            /* Go on to remove the redundant subspace from global pools after stopping this thread. */
            void ThreadStopDelete_Pool(int& target_depth);
            bool compare_sequence(vector<int>& sequence, int& target_depth);
    public:
        /* Takes config information and defines all runtime parameters from those strings. */
        void assign_parameter(vector<string> setting);
        /* Initial function that does all the necessary sequential preprocessing. */
        void solve(string filename,int thread_num);
};

#endif
