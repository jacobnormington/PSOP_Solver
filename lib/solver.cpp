#include "solver.hpp"
#include "hash.hpp"
#include "precedence.hpp"
extern "C" {
#include "LKH/LKHmain.h"
}
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <cmath>
#include <chrono>
#include <math.h>
#include <limits>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>
#include <vector>
#include <thread>
#include <queue>
#include <deque>
#include <mutex>
#include <bits/stdc++.h>
#include <unordered_map>
#include <condition_variable>

using namespace std;
#define TABLE_SIZE 541065431
#define TIME_FRAME 0.3
#define DATA_COLLECTION_START 5
#define RESTART_CNT_LIMIT 5
#define MAXIMUM_RECURCALL_CNT 50

//Restart Constants
static float exploitation_per = 0;
static int group_sample_time = -1;
static vector<int_64> selectT_arr;
static float inhis_mem_limit = -1;

int* bestBB_tour = NULL;
int best_cost = INT_MAX;
static Hash_Map history_table(TABLE_SIZE);
static vector<int> best_solution;
static Global_Pool GPQ;

//Variable for locking;
pthread_mutex_t Sol_lock = PTHREAD_MUTEX_INITIALIZER;
static mutex GPQ_lock, Split_lock;
static mutex asssign_mutex;
static mutex thread_load_mutex;
static condition_variable Idel;
static condition_variable Thread_Stop_Check;
static condition_variable Resume_State;
static vector<int> selected_orgin;
static mutex Select_Mutex;
static mutex Select_SharedMutex;
static mutex Resume_Lock;
static mutex launch_lck;

static atomic<int> selected_thread (-1);
static atomic<int> restart_cnt (0);
static atomic<int> total_restarts (0);
static atomic<int> active_thread (0);
static atomic<unsigned> idle_counter (0);
static atomic<size_t> resload_cnt (0);
static atomic<bool> time_out (false);
static atomic<bool> limit_insert (false);
static atomic<bool> check_status_safe (true);
static atomic<bool> resume_success (true);
static atomic<bool> resume_check (false);
static atomic<bool> exploit_init (false);

///////////Data Collection//////////////
static vector<int_64> num_resume;
static vector<int_64> num_stop;

///////////Thread Stop//////////////////
static deque<request_packet> request_buffer;
static mutex buffer_lock;
static mutex pause_lock;
static mutex ptselct_lock;
static atomic<int> stop_cnt (0);
static atomic<bool> stop_sig (false);
///////////////////////////////////////

///////////Thread Resume///////////////
static atomic<int> lowest_computedlb (INT_MAX);
static promise_stats restart_data;
///////////////////////////////////////

atomic<int> group_bestsolcnt(0);
atomic<int> refresh_threadcnt(0);

load_stats* thread_load;
bool end_heuristics = false;

//Config variable (Reading Only)
static string f_name;
static int global_pool_size = 0;
static int local_pool_size = 0;
static int tgroup_ratio = 0;
static int t_limit = 0;
static int thread_total = 0;
static int local_depth = 0;
static int max_edge_weight = 0;
static bool enable_workstealing = false;
static float pre_density = 0;
static bool enable_threadstop = false;
static bool enable_lkh = false;

//Shared resources
std::chrono::time_point<std::chrono::system_clock> start_time_limit;
std::chrono::time_point<std::chrono::system_clock> time_point;
static vector<vector<int>> dependent_graph;
static vector<vector<int>> outgoing_graph;
static vector<vector<edge>> in_degree;
static vector<vector<edge>> hung_graph;
static vector<vector<edge>> cost_graph;
static vector<Hungarian> initial_hungstate;
static sop_state default_state;
static vector<bool_64> abandon_wlk_array;
static vector<lptr_64> glp;
static vector<mutex_64> lp_lock;
static vector<unsigned_long_64> enumerated_nodes;

//Restart Variable
static int promise_Tlimit = 0;
static int promise_num = 0;
static int cut_off = 0;
static int res_id = 0;
static int group_cnt = 0;
static int grp_id = 0;
static int total_pcount = 0;
static vector<vector<Hungarian>> shared_lbstate;
static vector<condition_variable> restart_lbcompute;
static speed_resinfo** lbcompute_info = NULL;
static int_64* explored_cnt = NULL;
static int_64* explore_num = NULL;
static mutex_64* resinit_lock = NULL;
static mutex_64* group_lck = NULL;
static vector<atomwrapper<bool>> level_complete;
static vector<atomwrapper<bool>> level_processing;
static vector<atomwrapper<bool>> level_proceed;
static vector<vector<atomwrapper<bool>>> cfinish;
static atomic<int> psched_cnt(0);
static atomic<int> global_concentrate_lv(0);
static atomic<bool> lb_restart(false);
static vector<atomwrapper<bool>> busy_arr;

//LKH Variable
bool BB_Complete = false;
bool BB_SolFound = false;
bool local_searchinit = true;
bool initial_LKHRun = true;

//Test Variable
static vector<double> lp_time;
static vector<double> steal_wait;
static vector<vector<double>> proc_time;
static vector<int> steal_cnt;

thread LKH_thread;

void solver::assign_parameter(vector<string> setting) {
    t_limit = atoi(setting[0].c_str());
    cout << "Time limit = " << t_limit << endl;
    global_pool_size = atoi(setting[1].c_str());
    cout << "GPQ size = " << global_pool_size << endl;
    local_depth = atoi(setting[2].c_str());
    cout << "LPQ depth = " << local_depth << endl;

    inhis_mem_limit = atof(setting[3].c_str());
    cout << "History table mem limit = " << inhis_mem_limit << endl;

    exploitation_per = atof(setting[4].c_str())/float(100);
    cout << "Restart exploitation/exploration ratio is " << exploitation_per << endl;

    group_sample_time = atoi(setting[5].c_str());
    cout << "Group sample time = " << group_sample_time << endl;
    
    tgroup_ratio = atoi(setting[6].c_str());
    cout << "Number of promoising thread per exploitation group = " << tgroup_ratio << endl;

    if (!atoi(setting[7].c_str())) enable_workstealing = false;
    else enable_workstealing = true;

    if (!atoi(setting[8].c_str())) enable_threadstop = false;
    else enable_threadstop = true;

    if (!atoi(setting[9].c_str())) enable_lkh = false;
    else enable_lkh = true;

    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    problem_state.hungarian_solver.fix_row(src, dest);
	problem_state.hungarian_solver.fix_column(dest, src);
	problem_state.hungarian_solver.solve_dynamic();
    //number_of_lb++;
    return problem_state.hungarian_solver.get_matching_cost()/2;
}


bool solver::HistoryUtilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,HistoryNode** entry, int cost) {
    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket_location = (val + key.second) % TABLE_SIZE;
    HistoryNode* history_node = history_table.retrieve(key,bucket_location);

    if (history_node == NULL) return true;

    *found = true;
    *entry = history_node;
    HistoryContent content = history_node->Entry.load();
    *lowerbound = content.lower_bound;

    if (cost >= content.prefix_cost) return false;

    int target_ID = history_node->active_threadID;
    int imp = content.prefix_cost - cost;
    
    if (history_node->explored) {
        if (imp <= content.lower_bound - best_cost) {
            return false;
        }
        history_node->Entry.store({cost,content.lower_bound - imp});
        *lowerbound = content.lower_bound - imp;
        *entry = history_node;
        history_node->explored = false;
        history_node->active_threadID = thread_id;
    }
    else {
        if (enable_threadstop && active_thread > 0) {
            buffer_lock.lock();
            if (request_buffer.empty() || request_buffer.front().target_thread != target_ID || request_buffer.front().target_depth > (int)problem_state.cur_solution.size()) {
                //num_stop[thread_id].val++;
                request_buffer.push_front({problem_state.cur_solution.back(),(int)problem_state.cur_solution.size(),
                                           content.prefix_cost,bucket_location,target_ID});
                if (!stop_sig) {
                    stop_cnt = 0;
                    stop_sig = true;
                }
            }
            buffer_lock.unlock();
        }
        
        if (imp <= content.lower_bound - best_cost) {
            return false;
        }
        history_node->Entry.store({cost,content.lower_bound - imp});
        *lowerbound = content.lower_bound - imp;
        *entry = history_node;
        history_node->explored = false;
        history_node->active_threadID = thread_id;
    }
    
    return true;
}


bool Split_sort(const sop_state& src, const sop_state& dest) {
    if (src.cur_solution.size() == dest.cur_solution.size()) return src.load_info < dest.load_info;
    return src.cur_solution.size() < dest.cur_solution.size();
}

bool GPQ_sort(const instrct_node& src, const instrct_node& dest) {
    return src.load_info > dest.load_info;
}

bool Shared_pool_sort(const instrct_node& src, const instrct_node& dest) {
    if (src.sequence.size() == dest.sequence.size()) return src.load_info > dest.load_info;
    return src.sequence.size() < dest.sequence.size();
}

bool local_sort(const instrct_node& src, const instrct_node& dest) {
    return src.load_info > dest.load_info;
}

bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}

bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}

bool solver::Wlkload_Request() {
    bool terminate = true;

    if ((int)problem_state.cur_solution.size() == problem_state.initial_depth && active_thread > 0) {

        if (stop_init) {
            stop_init = false;
            stop_depth = -1;
            last_node = -1;
        }

        if (abandon_share) {
            abandon_share = false;
            ptselct_lock.lock();
            psched_cnt++;
            for (unsigned i = 0; i < selectT_arr.size(); i++) {
                if (thread_id == selectT_arr[i].val) {
                    selectT_arr[i].explored = true;
                    break;
                }
            }
            if (psched_cnt == total_pcount) {
                total_pcount = 0;
                psched_cnt = 0;
                cout << "Enabling Thread Restart" << endl;
                *GPQ.Promise.enable = true;
                selectT_arr.clear();
            }
            ptselct_lock.unlock();

            GPQ_lock.lock();
            if (!GPQ.Promise.is_empty()) sort(GPQ.Promise.space.back().begin(),GPQ.Promise.space.back().end(),Shared_pool_sort);
            cout << "Thread " << thread_id << " finished insertion, total promise node count is " << GPQ.Promise.space.back().size() << endl;
            GPQ_lock.unlock();
        }
        
        if (speed_search) {
            cout << "Speed restart disabled by thread " << thread_id << " upon finished with workload" << endl;
            speed_search = false;
            if (exploit_init) exploit_init = false;
        }

        if (concentrate_lv == 0) {

            thread_load_mutex.lock();
            thread_load[thread_id].out_of_work = true;
            thread_load[thread_id].data_cnt = 0;
            thread_load_mutex.unlock();

            if (abandon_work) initialize_node_sharing();

            terminate = Grab_from_GPQ(false);
            
            if (active_thread > 0 && terminate) {
                active_thread--;
                terminate = Steal_Workload();
                if (!terminate) active_thread++;
            }

            if (!terminate) {
                thread_load_mutex.lock();
                thread_load[thread_id].out_of_work = false;
                thread_load_mutex.unlock();
            }
        }
    }
    return terminate;
}

bool solver::Steal_Workload() {
    bool terminate = true;

    while (active_thread > 0 && wrksteal_pool.empty()) {
        idle_counter++;
        while (active_thread > 0 && wrksteal_pool.empty()) transfer_wlkload();
        idle_counter--;
    }

    if (!wrksteal_pool.empty()) {
        Generate_SolverState(wrksteal_pool.front());
        wrksteal_pool.pop_front();
        terminate = false;
    }

    return terminate;
}

void LB_Calculate(int g_id, int t_id) {
    cout << "Exploitation thread " << t_id << " in group " << g_id << " launched" << endl;

    while (exploit_init) {
        while (!cfinish[g_id][t_id].load() && exploit_init) continue;
        if (!exploit_init || time_out) {
            //cout << "Exploration thread " << t_id << " terminated" << endl;
            return;
        }

        cfinish[g_id][t_id].store(false);

        if (lbcompute_info[g_id][t_id].compute) {
            for (int i = lbcompute_info[g_id][t_id].start; i < lbcompute_info[g_id][t_id].finish; i++) {
                int w = (*lbcompute_info[g_id][t_id].lbproc_list)[i].n;
                
                shared_lbstate[g_id][t_id].fix_row(lbcompute_info[g_id][t_id].v, w);
                shared_lbstate[g_id][t_id].fix_column(w, lbcompute_info[g_id][t_id].v);
                shared_lbstate[g_id][t_id].solve_dynamic();
                (*lbcompute_info[g_id][t_id].lbproc_list)[i].lb = shared_lbstate[g_id][t_id].get_matching_cost()/2;

                shared_lbstate[g_id][t_id].undue_row(lbcompute_info[g_id][t_id].v,w);
                shared_lbstate[g_id][t_id].undue_column(w,lbcompute_info[g_id][t_id].v);
            }
        }
        
        group_lck[g_id].lck.lock();
        explored_cnt[g_id].val++;
        //if (lbcompute_info[g_id][t_id].compute) cout << "thread " << t_id << " finished computing LB explored cnt is " << explored_cnt[g_id].val << " total exploitation thread count is " << explore_num[g_id].val << endl;
        //else cout << "thread " << t_id << " did not explored " << endl;

        if (explored_cnt[g_id].val == explore_num[g_id].val) {
            //cout << "Enter Unlocking by thread " << t_id << endl;
            level_processing[g_id].store(true);
            while (!level_complete[g_id].load()) {
                //cout << "gid is " << g_id << endl;
                std::unique_lock<std::mutex> explore_guard(resinit_lock[g_id].lck);
                explore_guard.unlock();
                restart_lbcompute[g_id].notify_all();
            }
            level_proceed[g_id].store(true);
            //cout << "Unlocked by thread " << t_id << endl;
        }
        group_lck[g_id].lck.unlock();
    }
    
    return;
}

void solver::transfer_wlkload() {
    int designated_thread = -1;
    while (active_thread > 0 && wrksteal_pool.empty()) {
        designated_thread = rand() % thread_total;
        //cout << "designated thread is " << designated_thread << endl;
        if (busy_arr[designated_thread].load()) continue;
        lp_lock[designated_thread].lck.lock();
        busy_arr[thread_id].store(true);
        size_t transfer_size = ceil((1 - pre_density) * glp[designated_thread].local_pool->size());
        while (transfer_size > 0 && !glp[designated_thread].local_pool->empty()) {
            if (glp[designated_thread].local_pool->back().load_info < best_cost) {
                wrksteal_pool.push_back(glp[designated_thread].local_pool->back());
                transfer_size--;
            }
            *(glp[designated_thread].local_pool->back().invalid_ptr) = true;
            glp[designated_thread].local_pool->pop_back();
        }
        busy_arr[thread_id].store(false);
        lp_lock[designated_thread].lck.unlock();
    }
    return;
}

void solver::StopCurThread(int target_depth) {
    if (!stop_init) {
        stop_init = true;
        stop_depth = target_depth;
        ThreadStopDelete_Pool(target_depth);
    }
    else if (target_depth < stop_depth) {
        stop_depth = target_depth;
        ThreadStopDelete_Pool(target_depth);
    }
    return;
}

void solver::CheckStop_Request() {

    if (!request_buffer.empty() && !thread_load[thread_id].stop_checked) {
        buffer_lock.lock();
        if (request_buffer.empty()) {
            stop_sig = false;
            stop_cnt = 0;
            for (int i = 0; i < thread_total; i++) thread_load[i].stop_checked = false;
            buffer_lock.unlock();
            return;
        }
        int target_prefix_cost = request_buffer.back().target_prefix_cost;
        int target_depth = request_buffer.back().target_depth;
        int target_lastnode = request_buffer.back().target_lastnode;
        int target_key = request_buffer.back().target_key;
        buffer_lock.unlock();
        
        if (thread_stop_check(target_prefix_cost, target_depth, target_lastnode, target_key)) StopCurThread(target_depth);

        thread_load[thread_id].stop_checked = true;
        stop_cnt++;

        pause_lock.lock();
        if (stop_cnt >= active_thread) {
            buffer_lock.lock();
            request_buffer.pop_back();
            for (int i = 0; i < thread_total; i++) thread_load[i].stop_checked = false;
            stop_cnt = 0;
            if (request_buffer.empty()) stop_sig = false;
            else stop_sig = true;
            buffer_lock.unlock();
        }
        pause_lock.unlock();
    }

    return;
}

void solver::ThreadStopDelete_Pool(int& target_depth) {

    GPQ_lock.lock();

    for (unsigned k = 0; k < GPQ.Abandoned.size(); k++) {
        if (!GPQ.Abandoned[k].deprecated && compare_sequence(GPQ.Abandoned[k].sequence,target_depth) && GPQ.Abandoned[k].partial_active_path.check_deprecation_status(target_depth - 1)) {
            GPQ.Abandoned[k].partial_active_path.incre_children_cnt(Allocator);
            if (GPQ.Abandoned[k].root_his_node != NULL) {
                GPQ.Abandoned[k].root_his_node->explored = true;
            }
            GPQ.Abandoned[k].deprecated = true;
        }
    }

    for (unsigned i = 0; i < GPQ.Promise.space.size(); i++) {
        for (unsigned k = 0; k < GPQ.Promise.space[i].size(); k++) {
            if (!GPQ.Promise.space[i][k].deprecated && compare_sequence(GPQ.Promise.space[i][k].sequence,target_depth) && GPQ.Promise.space[i][k].partial_active_path.check_deprecation_status(target_depth - 1)) {
                GPQ.Promise.space[i][k].partial_active_path.incre_children_cnt(Allocator);
                if (GPQ.Promise.space[i][k].root_his_node != NULL) {
                    GPQ.Promise.space[i][k].root_his_node->explored = true;
                }
                GPQ.Promise.space[i][k].deprecated = true;
            }
        }
    }

    GPQ_lock.unlock();

    return;
}

void solver::regenerate_hungstate() {
    problem_state.hungarian_solver.n = initial_hungstate[thread_id].n;
    problem_state.hungarian_solver.N = initial_hungstate[thread_id].N;
    problem_state.hungarian_solver.max_edge_weight = initial_hungstate[thread_id].max_edge_weight;
    problem_state.hungarian_solver.max_match = initial_hungstate[thread_id].max_match;
    problem_state.hungarian_solver.cost = initial_hungstate[thread_id].cost;
    problem_state.hungarian_solver.lx = initial_hungstate[thread_id].lx;
    problem_state.hungarian_solver.ly = initial_hungstate[thread_id].ly;
    problem_state.hungarian_solver.xy = initial_hungstate[thread_id].xy;
    problem_state.hungarian_solver.yx = initial_hungstate[thread_id].yx;
    problem_state.hungarian_solver.S = initial_hungstate[thread_id].S;
    problem_state.hungarian_solver.T = initial_hungstate[thread_id].T;
    problem_state.hungarian_solver.slack = initial_hungstate[thread_id].slack;
    problem_state.hungarian_solver.slackx = initial_hungstate[thread_id].slackx;
    problem_state.hungarian_solver.prev = initial_hungstate[thread_id].prev;
    problem_state.hungarian_solver.fixed_rows = initial_hungstate[thread_id].fixed_rows;
    problem_state.hungarian_solver.fixed_columns = initial_hungstate[thread_id].fixed_columns;
    return;
}

void solver::Generate_SolverState(instrct_node& sequence_node) {
    int taken_node = -1;
    int counter = 0;
    int size = sequence_node.sequence.size();
    int cur_node = sequence_node.sequence.front();

    problem_state.cur_solution.clear();
    problem_state.cur_cost = 0;
    problem_state.taken_arr = default_state.taken_arr;
    problem_state.depCnt = default_state.depCnt;

    problem_state.cur_solution.push_back(0);
    problem_state.taken_arr[cur_node] = 1;
    for (int vertex : dependent_graph[cur_node]) problem_state.depCnt[vertex]--;
    counter++;
    
    regenerate_hungstate();

    while (counter < size) {
        taken_node = sequence_node.sequence[counter];
        problem_state.taken_arr[taken_node] = 1;
        for (int vertex : dependent_graph[taken_node]) problem_state.depCnt[vertex]--;
        problem_state.cur_cost += cost_graph[cur_node][taken_node].weight;
        problem_state.cur_solution.push_back(taken_node);
        problem_state.hungarian_solver.fix_row(cur_node, taken_node);
        problem_state.hungarian_solver.fix_column(taken_node, cur_node);
        cur_node = taken_node;
        counter++;
    }

    problem_state.initial_depth = problem_state.cur_solution.size();
    problem_state.load_info = sequence_node.load_info;
    problem_state.originate = sequence_node.originate;
    current_hisnode = sequence_node.root_his_node;
    
    cur_active_tree.generate_path(sequence_node.partial_active_path);
    cur_active_tree.set_threadID(thread_id, thread_total);
    if (current_hisnode != NULL && !current_hisnode->explored) {
        current_hisnode->active_threadID = thread_id;
    }

    boost::dynamic_bitset<> bit_vector(node_count, false);
    for (auto node : problem_state.cur_solution) {
        bit_vector[node] = true;
    }
    int last_element = problem_state.cur_solution.back();
    problem_state.key = make_pair(bit_vector,last_element);

    lb_curlv = sequence_node.load_info;

    return;
}

bool solver::assign_workload(node& transfer_node, pool_dest destination, space_ranking problem_property, HistoryNode* temp_hisnode) {
    int taken_n = transfer_node.n;
    int lb = transfer_node.lb;

    if (lb >= best_cost) {
        cur_active_tree.incre_children_cnt(Allocator);
        if (temp_hisnode != NULL && temp_hisnode->active_threadID == thread_id) {
            temp_hisnode->explored = true;
        }
        return false;
    }

    vector<int> target_solution = problem_state.cur_solution;
    target_solution.push_back(taken_n);
    instrct_node workload = instrct_node(target_solution,problem_state.originate+1,lb,-1,
                                         cur_active_tree,temp_hisnode);

    switch (destination) {
        case GLOBAL_POOL:
            switch (problem_property) {
                case ABANDONED:
                    GPQ_lock.lock();
                    GPQ.Abandoned.push_back(workload);
                    GPQ_lock.unlock();
                    break;
                case PROMISING:
                    GPQ_lock.lock();
                    GPQ.Promise.space.back().push_back(workload);                                                                                         
                    GPQ_lock.unlock();
                    break;
                case UNEXPLORED:
                    GPQ_lock.lock();
                    GPQ.Unknown.push_back(workload);                                     
                    GPQ_lock.unlock();
                    break;
            }
            break;
        case LOCAL_POOL:
            workload.invalid_ptr = &(transfer_node.node_invalid);
            local_pool->push_back(workload);
            break;
    }
    return true;
}

bool solver::push_to_pool(pool_dest decision, space_ranking problem_property) {
    bool insertion = false;

    lp_lock[thread_id].lck.lock();
    while (!local_pool->empty()) {
        if (local_pool->back().load_info >= best_cost) {
            local_pool->back().partial_active_path.incre_children_cnt(Allocator);
            if (local_pool->back().root_his_node != NULL && local_pool->back().root_his_node->active_threadID == thread_id) {
                local_pool->back().root_his_node->explored = true;
            }
            *(local_pool->back().invalid_ptr) = true;
            local_pool->pop_back();
        }
        else {
            switch(decision) {
                case GLOBAL_POOL:
                    switch (problem_property) {
                        case ABANDONED:
                            GPQ_lock.lock();
                            GPQ.Abandoned.push_back(local_pool->back());
                            GPQ_lock.unlock();
                            break;
                        case PROMISING:
                            GPQ_lock.lock();
                            GPQ.Promise.space.back().push_front(local_pool->back());
                            GPQ_lock.unlock();
                            break;
                        case UNEXPLORED:
                            GPQ_lock.lock();
                            GPQ.Unknown.push_back(local_pool->back());
                            GPQ_lock.unlock();
                            break;
                    }
                    break;
                case LOCAL_POOL:
                    break;
            }
            *(local_pool->back().invalid_ptr) = true;
            local_pool->pop_back();
            insertion = true;
        }
    }
    lp_lock[thread_id].lck.unlock();

    return insertion;
}


void solver::notify_finished() {
    std::unique_lock<std::mutex> idel_lck(Split_lock);
    idel_lck.unlock();
    Idel.notify_all();
    return;
}

void solver::push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,HistoryNode** entry,bool backtracked) {
    if (entry == NULL) history_table.insert(key,problem_state.cur_cost,lower_bound,thread_id,backtracked,problem_state.cur_solution.size());
    else *entry = history_table.insert(key,problem_state.cur_cost,lower_bound,thread_id,backtracked,problem_state.cur_solution.size());
    return;
}

bool solver::Split_level_check(deque<sop_state>* solver_container) {
    unsigned target_level = solver_container->front().cur_solution.size();
    for (auto node : *solver_container) {
        if (node.cur_solution.size() != target_level) return true;
    }

    return false;
}

bool solver::Is_promisethread(int t_id) {
    if (selectT_arr.empty()) return false;

    ptselct_lock.lock();
    for (int_64 id_val : selectT_arr) {
        if (id_val.val == t_id && id_val.explored == false) {
            ptselct_lock.unlock();
            return true;
        }
    }
    ptselct_lock.unlock();
    return false;
}

void solver::Select_GroupMember(int unpromise_cnt) {

    if (!total_pcount && GPQ.Promise.is_empty()) {
        for (int k = 0; k < thread_total; k++) {
            if (GPQ.Unknown.empty()) {
                thread_load[k].next_action = abandon_action::CONTINUE_EXPLORATION;
                abandon_wlk_array[k].val = false;
                cout << "Thread " << k  << " continues exploration" << endl;
            }
            else {
                if (thread_load[k].space_ranking == space_state::NOT_PROMISING) {
                    thread_load[k].next_action = abandon_action::EXPLORE_UNKNOWN;
                    abandon_wlk_array[k].val = true;
                    cout << "Thread " << k << " explore unknown" << endl;
                }
                else {
                    thread_load[k].next_action = abandon_action::CONTINUE_EXPLORATION;
                    abandon_wlk_array[k].val = false;
                    cout << "Thread " << k  << " with gb_cnt == " << thread_load[k].data_cnt << " continues exploration" << endl;
                }
            }
        }
    }
    else {  
        int exploitation_cnt = ceil((float)unpromise_cnt*exploitation_per);

        if (!group_bestsolcnt) {
            for (int i = exploitation_cnt; i > 0; i--) {
                if (i % selectT_arr.size() == 0) {
                    exploitation_cnt = i;
                    break;
                }
            }
        }

        int exploration_cnt = unpromise_cnt - exploitation_cnt;

        if (!group_bestsolcnt) {
            res_id = 0;
            group_cnt = 0;
            grp_id = 0;
            cut_off = exploitation_cnt / selectT_arr.size();
            exploit_init = true;
            shared_lbstate = vector<vector<Hungarian>>(promise_num);
            for (int i = 0; i < promise_num; i++) shared_lbstate[i] = vector<Hungarian>(cut_off);
            lbcompute_info = new speed_resinfo*[promise_num];
            for (int i = 0; i < promise_num; i++) {
                lbcompute_info[i] = new speed_resinfo[cut_off];
            }
            restart_lbcompute = vector<condition_variable>(promise_num);
            resinit_lock = new mutex_64[promise_num];
            explore_num = new int_64[promise_num];
            for (int i = 0; i < promise_num; i++) explore_num[i].val = cut_off;
            explored_cnt = new int_64[promise_num];
            for (int i = 0; i < promise_num; i++) explored_cnt[i].val = 0;
            group_lck = new mutex_64[promise_num];
            cfinish = vector<vector<atomwrapper<bool>>>(promise_num);
            for (int i = 0; i < promise_num; i++) {
                cfinish[i] = vector<atomwrapper<bool>>(cut_off);
            }
            level_complete = vector<atomwrapper<bool>>(promise_num);
            level_processing = vector<atomwrapper<bool>>(promise_num);
            level_proceed = vector<atomwrapper<bool>>(promise_num);
            cout << exploitation_cnt << " threads are selected for beam restart" << endl;
        }

        for (int k = 0; k < thread_total; k++) {
            if (thread_load[k].space_ranking == space_state::NOT_PROMISING) {
                if (exploitation_cnt > 0) {
                    thread_load[k].next_action = abandon_action::EXPLOITATION_PROMISE;
                    if (!Is_promisethread(k)) {
                        abandon_wlk_array[k].val = true;
                        exploitation_cnt--;
                        cout << "Thread " << k << " explore promising" << endl;
                    }
                }
                else if (exploration_cnt > 0) {
                    thread_load[k].next_action = abandon_action::EXPLORE_UNKNOWN;
                    if (group_bestsolcnt && !Is_promisethread(k)) {
                        abandon_wlk_array[k].val = true;
                        exploration_cnt--;
                        cout << "Thread " << k << " explore unknown" << endl;
                    }
                    else {
                        abandon_wlk_array[k].val = false;
                        cout << "Thread " << k << " continues exploration" << endl;
                    }
                }
            }
            else {
                thread_load[k].next_action = abandon_action::CONTINUE_EXPLORATION;
                cout << "Thread " << k  << " with gb_cnt == " << thread_load[k].data_cnt << " keeps exploring" << endl;
            }
        }
    }
    return;
}

void solver::Thread_Selection() {
    int max = 0;
    int id = -1;
    int unpromise_cnt = 0;
    vector<int> restart_info;
    vector<int_64> temp_select;
    restart_info = vector<int>(thread_total, -1);

    thread_load_mutex.lock();

    ptselct_lock.lock();
    selectT_arr.clear();
    ptselct_lock.unlock();

    for (int i = 0; i < thread_total; i++) {
        if (!thread_load[i].data_cnt) {
            thread_load[i].space_ranking = space_state::NOT_PROMISING;
            unpromise_cnt++;
        }
        else if (thread_load[i].data_cnt >= 1) thread_load[i].space_ranking = space_state::IS_PROMISING;
    }

    if (!unpromise_cnt) {
        cout << "No unpromise threads. Search progress is good skipping to next cycle" << endl;
        thread_load_mutex.unlock();
        return;
    }

    for (int i = 0; i < thread_total; i++) {
        restart_info[i] = thread_load[i].data_cnt;
        cout << restart_info[i] << ",";
    }
    cout << endl;

    cout << "Selected threads are --> ";

    double src_found_time = 0;
    double dest_found_time = 0;

    for (int i = 0; i < promise_Tlimit; i++) {
        for (int j = 0; j < thread_total; j++) {
            if (restart_info[j] > 1) {
                if (restart_info[j] > max) {
                    max = restart_info[j];
                    id = j;
                }
                else if (restart_info[j] == max) {
                    src_found_time = std::chrono::duration<double>(thread_load[id].found_time - start_time_limit).count();
                    dest_found_time = std::chrono::duration<double>(thread_load[j].found_time - start_time_limit).count();
                    if (dest_found_time > src_found_time) id = j;
                }
            }
        }
        if (id != -1) {
            int_64 target;
            target.val = id;
            target.explored = false;
            restart_info[id] = -1;
            cout << target.val << ",";
            temp_select.push_back(target);
        }
        max = 0;
        id = -1;
    }
    cout << endl;

    
    if (!temp_select.empty() && group_bestsolcnt) {
        GPQ_lock.lock();
        GPQ.Promise.space.push_back(deque<instrct_node>());
        global_concentrate_lv++;
        cout << "current promise pool count is " << global_concentrate_lv << endl;
        GPQ_lock.unlock();
    }

    total_pcount = temp_select.size();

    ptselct_lock.lock();
    for (unsigned i = 0; i < temp_select.size(); i++) {
        selectT_arr.push_back(temp_select[i]);
    }
    ptselct_lock.unlock();
    
    promise_num = selectT_arr.size();
    Select_GroupMember(unpromise_cnt);

    for (int i = 0; i < thread_total; i++) {
        thread_load[i].data_cnt = 0;
    }
    thread_load_mutex.unlock();
    return;
}


void solver::initialize_node_sharing() {

    if (((thread_load[thread_id].next_action == abandon_action::EXPLOITATION_PROMISE && *GPQ.Promise.enable) ||
        (thread_load[thread_id].next_action == abandon_action::EXPLORE_UNKNOWN)) && abandon_work) {
        abandon_work = false;
       
        thread_load_mutex.lock();
        int backup_gb_cnt = thread_load[thread_id].data_cnt;
        thread_load_mutex.unlock();
        
        recur_stack.push_back({problem_state,cur_active_tree});
        HistoryNode* backup_curhisnode = current_hisnode;
        bool Grabbed = false;
        while (true) {
            GPQ_lock.lock();
            if (    (thread_load[thread_id].next_action == abandon_action::EXPLORE_UNKNOWN && !GPQ.Unknown.empty())
                 || (thread_load[thread_id].next_action == abandon_action::EXPLOITATION_PROMISE && !GPQ.Promise.is_empty())) {
                Grabbed = grab_restartnode();
            }
            else {
                if (thread_load[thread_id].next_action == abandon_action::EXPLOITATION_PROMISE && GPQ.Promise.is_empty()) {
                    *GPQ.Promise.enable = false;
                    cout << "Thread " << thread_id << " exits sharing execution" << endl;
                }
                abandon_work = false;
                abandon_share = false;
                problem_state = recur_stack.back().state;
                cur_active_tree.generate_path(recur_stack.back().atree);
                recur_stack.pop_back();
                current_hisnode = backup_curhisnode;
                GPQ_lock.unlock();
                return;
            }
            GPQ_lock.unlock();

            if (Grabbed) {
                //Change variables for second restart
                
                thread_load_mutex.lock();
                thread_load[thread_id].data_cnt = 0;
                thread_load_mutex.unlock();

                concentrate_lv++;
                enumerate();
                concentrate_lv--;

                thread_load_mutex.lock();
                //thread_load[thread_id].data_cnt = backup_gb_cnt;
                thread_load[thread_id].data_cnt = backup_gb_cnt;
                thread_load_mutex.unlock();

                if (abandon_work) abandon_work = false;
            }
            auto cur_time = std::chrono::system_clock::now();

            if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
                time_out = true;
                active_thread = 0;
                problem_state = recur_stack.back().state;
                cur_active_tree.generate_path(recur_stack.back().atree);
                recur_stack.pop_back();
                current_hisnode = backup_curhisnode;
                return;
            }
            Grabbed = false;
        }
        problem_state = recur_stack.back().state;
        cur_active_tree.generate_path(recur_stack.back().atree);
        recur_stack.pop_back();
        current_hisnode = backup_curhisnode;
    }

    return;
}


void solver::Check_Restart_Status(deque<node>& enumeration_list, deque<node>& curlocal_nodes) {

    if (thread_id == 0) {
        auto elasped_time = std::chrono::duration<double>(std::chrono::system_clock::now() - time_point).count();
        
        if (!lb_restart && !group_bestsolcnt && elasped_time > 2) {
            cout << "No solution detected. Issuing special LB restart" << endl;
            lb_restart = true;
	        Thread_Selection();
            time_point = std::chrono::system_clock::now();
        }
        else if (group_bestsolcnt && psched_cnt == total_pcount && elasped_time > group_sample_time) {
            abandon_share = false;
            abandon_work = false;
            cout << "examined by thread 0" << endl;
            Thread_Selection();
            time_point = std::chrono::system_clock::now();
        }
    }
    
    if (!speed_search && Is_promisethread(thread_id)) {
        abandon_share = true;
    }
    
    if (abandon_wlk_array[thread_id].val) {
        //cout << "Thread " << thread_id << " obtained sharing status" << endl;
        abandon_work = true;
        abandon_wlk_array[thread_id].val = false;
    }

    if (concentrate_lv > MAXIMUM_RECURCALL_CNT) {
        if (abandon_share) {
            ptselct_lock.lock();
            psched_cnt++;
            for (unsigned i = 0; i < selectT_arr.size(); i++) {
                if (thread_id == selectT_arr[i].val) {
                    selectT_arr[i].explored = true;
                    break;
                }
            }
            if (psched_cnt == total_pcount) {
                psched_cnt = 0;
                cout << "Enabling Thread Restart" << endl;
                GPQ_lock.lock();
                if (!GPQ.Promise.is_empty()) *GPQ.Promise.enable = true;
                GPQ_lock.unlock();
                selectT_arr.clear();
            }
            ptselct_lock.unlock();
            abandon_share = false;
        }
        if (abandon_work) abandon_work = false;
    }
    
    if (abandon_share && !abandon_work) {
        if (!group_bestsolcnt) {
            cout << "Beam restart initialized on thread " << thread_id << endl;
            mg_id = get_mgid();
            //cout << "master group id of promising thread " << thread_id << " is " << mg_id << endl;
            for (unsigned i = 0; i < shared_lbstate[mg_id].size(); i++) {
                shared_lbstate[mg_id][i] = problem_state.hungarian_solver;
            }
            speed_search = true;
            abandon_share = false;
            ptselct_lock.lock();
            psched_cnt++;
            for (unsigned i = 0; i < selectT_arr.size(); i++) {
                if (thread_id == selectT_arr[i].val) {
                    selectT_arr[i].explored = true;
                    break;
                }
            }
            if (psched_cnt == total_pcount) {
                total_pcount = 0;
                psched_cnt = 0;
                cout << "Enabling Thread Restart" << endl;
                selectT_arr.clear();
            }
            ptselct_lock.unlock();
        }
        else {
            push_to_pool(pool_dest::GLOBAL_POOL,space_ranking::PROMISING);
            while (!enumeration_list.empty()) {
                assign_workload(enumeration_list.back(),pool_dest::GLOBAL_POOL,space_ranking::PROMISING,enumeration_list.back().his_entry);
                enumeration_list.pop_back();
            }
        }
    }
    else if (!abandon_share && abandon_work) {
        if (!group_bestsolcnt) {
            int t_id = -1;
            int g_id = -1;
            abandon_work = false;
            cout << "Thread " << thread_id << " launches lb sharing" << endl;
            launch_lck.lock();
            t_id = res_id;
            g_id = grp_id;
            res_id++;
            group_cnt++;
            if (group_cnt == cut_off) {
                group_cnt = 0;
                res_id = 0;
                grp_id++;
            }
            launch_lck.unlock();
            LB_Calculate(g_id, t_id);
        }
        else {
            push_to_pool(pool_dest::GLOBAL_POOL,space_ranking::ABANDONED);
            while (!enumeration_list.empty()) {
                assign_workload(enumeration_list.back(),pool_dest::GLOBAL_POOL,space_ranking::ABANDONED,enumeration_list.back().his_entry);
                enumeration_list.pop_back();
            }
        }
    }
    return;
}

bool solver::Check_Local_Pool(deque<node>& enumeration_list,deque<node>& curlocal_nodes) {
    bool pushed_to_local = false;
    bool insert_condition = false;

    lp_lock[thread_id].lck.lock();

    if (idle_counter <= 0) insert_condition = local_pool->size() < 0.2 * (size_t)local_pool_size;
    else insert_condition = local_pool->size() <  0.8 * (size_t)local_pool_size;

    if (insert_condition && !enumeration_list.empty()) {
        if (problem_state.cur_solution.size() <= float(local_depth) / 100 * node_count || idle_counter > 0) {
            while (!enumeration_list.empty()) {
                curlocal_nodes.push_front(enumeration_list.back());
                enumeration_list.pop_back();
                if (!assign_workload(curlocal_nodes.front(),pool_dest::LOCAL_POOL,space_ranking::UNEXPLORED,curlocal_nodes.front().his_entry)) {
                    curlocal_nodes.pop_front();
                }
            }
        }
        if (!local_pool->empty()) {
            sort(local_pool->begin(),local_pool->end(),local_sort);
        }
    }

    lp_lock[thread_id].lck.unlock();

    return pushed_to_local;
}

bool solver::EnumerationList_PreProcess(deque<node>& enumeration_list,deque<node>& curlocal_nodes) {
    if (enumeration_list.back().lb >= best_cost 
        || stop_init
        || (enable_threadstop && enumeration_list.back().his_entry != NULL 
                              && enumeration_list.back().his_entry->Entry.load().prefix_cost < enumeration_list.back().partial_cost)
       )
    {   
        cur_active_tree.incre_children_cnt(Allocator);
        if (enumeration_list.back().his_entry != NULL && enumeration_list.back().his_entry->active_threadID == thread_id) {
            enumeration_list.back().his_entry->explored = true;
        }
        enumeration_list.pop_back();
        //current_stored_node--;
        
        if (enumeration_list.empty() && !curlocal_nodes.empty()) retrieve_from_local(curlocal_nodes,enumeration_list);

        return true;
    }
    return false;
}


bool solver::thread_stop_check(int target_prefix_cost, int target_depth, int target_lastnode, int target_key) {

    if ((int)problem_state.cur_solution.size() < target_depth) return false;
    else if (problem_state.cur_solution[target_depth-1] != target_lastnode) return false;

    if (cur_active_tree.get_element(target_depth - 1)->deprecated) return true;

    int taken_node = -1;
    int cur_node = problem_state.cur_solution.front();
    int computed_cost = 0;
    boost::dynamic_bitset<> bit_vector(node_count);
    for (int k = 0; k < target_depth; k++) {
        taken_node = problem_state.cur_solution[k];
        computed_cost += cost_graph[cur_node][taken_node].weight;
        bit_vector[problem_state.cur_solution[k]] = true;
        cur_node = taken_node;
    }
    auto key = make_pair(bit_vector,target_lastnode);
    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket_location = (val + key.second) % TABLE_SIZE;

    if (bucket_location != target_key) return false;
    else if (computed_cost <= target_prefix_cost) return false;

    cur_active_tree.get_element(target_depth - 1)->deprecated = true;

    return true;
}

bool solver::compare_sequence(vector<int>& sequence, int& target_depth) {
    if (sequence.size() < 1) return false;
    if (sequence.size() <= (unsigned)target_depth) return false;
    return true;
}

int solver::get_mgid() {
    ptselct_lock.lock();
    for (unsigned i = 0; i < selectT_arr.size(); i++) {
        if (selectT_arr[i].val == thread_id) {
            ptselct_lock.unlock();
            return i;
        }
    }
    ptselct_lock.unlock();
    return -1;
}

void solver::enumerate() {
    if (time_out) return;

    while (true) {
        deque<node> ready_list;
        deque<node> curlocal_nodes;
        
        for (int i = node_count-1; i >= 0; i--) {
            if (!problem_state.depCnt[i] && !problem_state.taken_arr[i]) {
                //Push vertices with 0 depCnt into the ready list
                ready_list.push_back(node(i,-2));
            }
        }

        int last_element = problem_state.cur_solution.back();
        int taken_node = 0;
        int u = 0;
        int v = 0;

        //cout << "Lv " << problem_state.cur_solution.size() << " ";
        
        deque<node> enumeration_list;
        deque<node> lbproc_list;
        bool limit_insertion = false;

        for (int i = 0; i < (int)ready_list.size(); i++) {
            node dest = ready_list[i];
            int src = problem_state.cur_solution.back();
            problem_state.cur_solution.push_back(dest.n);
            problem_state.cur_cost += cost_graph[src][dest.n].weight;
            if (problem_state.cur_cost >= best_cost) {
                problem_state.cur_solution.pop_back();
                problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                continue;
            }
            if (problem_state.cur_solution.size() == (size_t)node_count) {
                if (problem_state.cur_cost < best_cost) {
                    pthread_mutex_lock(&Sol_lock);
                    if (problem_state.cur_cost < best_cost) {
                        best_solution = problem_state.cur_solution;
                        best_cost = problem_state.cur_cost;
                        if (enable_lkh) {
                            for (int i = 0; i < node_count; i++) bestBB_tour[i] = best_solution[i] + 1;
                            BB_SolFound = true;
                        }

                        if (!group_bestsolcnt) {
                            thread_load_mutex.lock();
                            for (int i = 0; i < thread_total; i++) thread_load[i].data_cnt = 0;
                            thread_load_mutex.unlock();
                        }
                        group_bestsolcnt++;

                        if (std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() > DATA_COLLECTION_START) {
                            thread_load_mutex.lock();
                            thread_load[thread_id].data_cnt++;
                            thread_load[thread_id].found_time = std::chrono::system_clock::now();
                            thread_load_mutex.unlock();
                        }
                        
                        cout << "Best Cost = " << best_cost << " Found in Thread " << thread_id;
                        cout << " at time = " << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << " [" << thread_load[thread_id].data_cnt << " ]" << endl;
                        
                        /*
                        cout << "current solution is ";
                        for (auto node : problem_state.cur_solution) cout << node << ",";
                        cout <<  endl;
                        */
                        
                        if (speed_search) {
                            cout << "Speed restart disabled by thread " << thread_id << endl;
                            speed_search = false;
                            exploit_init = false;
                        }
                    }        
                    pthread_mutex_unlock(&Sol_lock);
                }
                problem_state.suffix_cost = problem_state.cur_cost;
                problem_state.cur_solution.pop_back();
                problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                continue;
            }
            problem_state.cur_solution.pop_back();
            problem_state.cur_cost -= cost_graph[src][dest.n].weight;
            lbproc_list.push_back(dest);
        }

        //auto start_time_wait = std::chrono::system_clock::now();

        //Issue LB sharing
        while (speed_search) {
            bool complete = false;
            bool under_load = false;
            int total_load = lbproc_list.size();
            int load_size = total_load / shared_lbstate[mg_id].size();
            int loc = 0;
            int load_cnt = total_load;

            if (!exploit_init) {
                cout << "Speed restart disabled by thread " << thread_id << endl;
                speed_search = false;
                exploit_init = false;
                break;
            }
            
            if (load_size == 0 && total_load < (int)shared_lbstate[mg_id].size()) {
                cout << "Level " << problem_state.cur_solution.size() << " underload" << endl;
                under_load = true;
            }

            for (int i = 0; i < cut_off; i++) {
                if (loc + load_size > total_load) load_size = total_load - loc;
                lbcompute_info[mg_id][i].u = problem_state.cur_solution.end()[-2];
                lbcompute_info[mg_id][i].v = problem_state.cur_solution.back();
                lbcompute_info[mg_id][i].start = loc;
                lbcompute_info[mg_id][i].finish = loc + load_size;
                lbcompute_info[mg_id][i].lbproc_list = &lbproc_list;
                load_cnt--;
                if (under_load && load_cnt < 0) lbcompute_info[mg_id][i].compute = false;
                else lbcompute_info[mg_id][i].compute = true;
                loc += load_size;
            }
            
            
            for (int i = 0; i < cut_off; i++) cfinish[mg_id][i].store(true);
            level_complete[mg_id].store(false);
            level_processing[mg_id].store(false);
            level_proceed[mg_id].store(false);
            if (problem_state.cur_solution.size() % 100 == 0) cout << "Level " << problem_state.cur_solution.size() << " initialized in group " << mg_id << endl;
            while (!level_processing[mg_id].load()) {
                std::unique_lock<std::mutex> explore_guard(resinit_lock[mg_id].lck);
                restart_lbcompute[mg_id].wait(explore_guard);
            }
            level_complete[mg_id].store(true);
            while (!level_proceed[mg_id].load()) continue;
            if (problem_state.cur_solution.size() % 100 == 0) cout << "Level " << problem_state.cur_solution.size() << " completed in group " << mg_id << endl;
            if (explored_cnt[mg_id].val == explore_num[mg_id].val) explored_cnt[mg_id].val = 0;
            else {
                cout << "error: not all threads have finished with the current workload " << endl;
                cout << "explored cnt == " << explored_cnt[mg_id].val << ", explore_num == " << explore_num[mg_id].val << endl;
                exit(1);
            }
            complete = true;
            if (complete) break;
        }

        for (int i = 0; i < (int)lbproc_list.size(); i++) {
            node dest = lbproc_list[i];
            int src = problem_state.cur_solution.back();
            problem_state.cur_solution.push_back(dest.n);
            problem_state.cur_cost += cost_graph[src][dest.n].weight;
            int lower_bound = -1;
            bool taken = false;
            problem_state.key.first[dest.n] = true;
            problem_state.key.second = dest.n;
            HistoryNode* his_node = NULL;
            Active_Node* active_node = NULL;

            if (speed_search && lbproc_list[i].lb > 0) {
                lower_bound = lbproc_list[i].lb;
                if (lbproc_list[i].lb >= best_cost) {
                    problem_state.cur_solution.pop_back();
                    problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                    problem_state.key.first[dest.n] = false;
                    problem_state.key.second = last_element;
                    continue;
                }
            }
            else {
                bool decision = HistoryUtilization(problem_state.key,&lower_bound,&taken,&his_node,problem_state.cur_cost);
                //enumerated_nodes[thread_id][problem_state.cur_solution.size()]++;

                if (!taken) {
                    lower_bound = dynamic_hungarian(src,dest.n);
                    if (history_table.get_cur_size() < inhis_mem_limit * history_table.get_max_size()) push_to_historytable(problem_state.key,lower_bound,&his_node,false);
                    else limit_insertion = true;
                    problem_state.hungarian_solver.undue_row(src,dest.n);
                    problem_state.hungarian_solver.undue_column(dest.n,src);
                }
                else if (taken && !decision) {
                    problem_state.cur_solution.pop_back();
                    problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                    problem_state.key.first[dest.n] = false;
                    problem_state.key.second = last_element;
                    continue;
                }

                if (lower_bound >= best_cost) {
                    if (his_node != NULL) {
                        HistoryContent content = his_node->Entry.load();
                        if (content.prefix_cost >= problem_state.cur_cost) his_node->explored = true;
                    }
                    problem_state.cur_solution.pop_back();
                    problem_state.cur_cost -= cost_graph[src][dest.n].weight;
                    problem_state.key.first[dest.n] = false;
                    problem_state.key.second = last_element;
                    continue;
                }
            }
            enumeration_list.push_back(lbproc_list[i]);
            enumeration_list.back().nc = cost_graph[src][dest.n].weight;
            enumeration_list.back().pushed = taken;
            enumeration_list.back().lb = lower_bound;
            enumeration_list.back().act_entry = active_node;
            enumeration_list.back().his_entry = his_node;
            enumeration_list.back().partial_cost = problem_state.cur_cost;
            problem_state.key.first[dest.n] = false;
            problem_state.key.second = last_element;

            problem_state.cur_solution.pop_back();
            problem_state.cur_cost -= cost_graph[src][dest.n].weight;
        }
        if (!enumeration_list.empty()) sort(enumeration_list.begin(),enumeration_list.end(),bound_sort);

        HistoryNode* history_entry = NULL;

        int lb_liminsert = lb_curlv;

        cur_active_tree.push_back(enumeration_list.size(),current_hisnode,Allocator);

        CheckStop_Request();

        while(!enumeration_list.empty()) {
            if(EnumerationList_PreProcess(enumeration_list,curlocal_nodes)) continue;

            //Check_Restart_Status(enumeration_list, curlocal_nodes);

            if (abandon_share || abandon_work) {
                curlocal_nodes.clear();
                break;
            }

            taken_node = enumeration_list.back().n;
            history_entry = enumeration_list.back().his_entry;
            lb_curlv = enumeration_list.back().lb;
            problem_state.key.first[taken_node] = true;
            problem_state.key.second = taken_node;
            
            if (group_bestsolcnt == 0 && problem_state.cur_solution.size() > restart_data.depth && enumeration_list.back().lb <= restart_data.lb) {
                restart_data.lb = enumeration_list.back().lb;
                restart_data.depth = problem_state.cur_solution.size();
                thread_load[thread_id].data_cnt++;
                thread_load[thread_id].found_time = std::chrono::system_clock::now();
            }

            enumeration_list.pop_back();

            if (!stop_init) {
                Check_Local_Pool(enumeration_list,curlocal_nodes);
            }
        
            u = problem_state.cur_solution.back();
            v = taken_node;
            problem_state.hungarian_solver.fix_row(u, v);
            problem_state.hungarian_solver.fix_column(v, u);
            problem_state.cur_cost += cost_graph[u][v].weight;

            if (speed_search) {
                for (unsigned i = 0; i < shared_lbstate[mg_id].size(); i++) {
                    shared_lbstate[mg_id][i].fix_row(u,v);
                    shared_lbstate[mg_id][i].fix_column(v,u);
                }
            }

            for (int vertex : dependent_graph[taken_node]) problem_state.depCnt[vertex]--;

            problem_state.cur_solution.push_back(taken_node);
            problem_state.taken_arr[taken_node] = 1;
            problem_state.suffix_cost = 0;

            HistoryNode* previous_hisnode = current_hisnode;
            current_hisnode = history_entry;

            enumerate();

            current_hisnode = previous_hisnode;

            for (int vertex : dependent_graph[taken_node]) problem_state.depCnt[vertex]++;
            problem_state.taken_arr[taken_node] = 0;
            problem_state.cur_solution.pop_back();
            problem_state.cur_cost -= cost_graph[u][v].weight;
            problem_state.hungarian_solver.undue_row(u,v);
            problem_state.hungarian_solver.undue_column(v,u);
            problem_state.key.first[taken_node] = false;
            problem_state.key.second = problem_state.cur_solution.back();

            if (speed_search) {
                for (unsigned i = 0; i < shared_lbstate[mg_id].size(); i++) {
                    shared_lbstate[mg_id][i].undue_row(u,v);
                    shared_lbstate[mg_id][i].undue_row(v,u);
                }
            }
            
            retrieve_from_local(curlocal_nodes,enumeration_list);
            
            if (thread_id == 0) {
                auto cur_time = std::chrono::system_clock::now();
                if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
                    time_out = true;
                    active_thread = 0;
                    return;
                }
            }
        }

        if (stop_init && (int)problem_state.cur_solution.size() <= stop_depth) {
            stop_init = false;
            stop_depth = -1;
            last_node = -1;
        }
        
        if (limit_insertion && history_table.get_cur_size() < history_table.get_max_size()) {
            push_to_historytable(problem_state.key,lb_liminsert,NULL,false);
        }

        cur_active_tree.pop_back(stop_init,Allocator);

        if (Wlkload_Request()) break;
    }

    return;
}

void solver::retrieve_from_local(deque<node>& curlocal_nodes, deque<node>& enumeration_list) {

    lp_lock[thread_id].lck.lock();
    if (!curlocal_nodes.empty()) {
        if (stop_init) {
            while (!curlocal_nodes.empty()) {
                if (curlocal_nodes.back().node_invalid) {
                    curlocal_nodes.pop_back();
                }
                else {
                    for (int k = local_pool->size() - 1; k >= 0; k--) {
                        if ((*local_pool)[k].invalid_ptr == &(curlocal_nodes.back().node_invalid)) {
                            //don't pop if it is still alive living in another thread node_grab
                            curlocal_nodes.pop_back();
                            (*local_pool)[k].partial_active_path.incre_children_cnt(Allocator);
                            if ((*local_pool)[k].root_his_node != NULL && (*local_pool)[k].root_his_node->active_threadID == thread_id) {
                                (*local_pool)[k].root_his_node->explored = true;
                            }
                            local_pool->erase(local_pool->begin() + k);
                        }
                    }
                }
            }
        }
        else {
            bool taken = false;
            while (!taken && !curlocal_nodes.empty()) {
                if (curlocal_nodes.back().node_invalid) {
                    curlocal_nodes.pop_back();
                }
                else {
                    taken = true;
                    enumeration_list.push_back(curlocal_nodes.back());
                    for (unsigned i = local_pool->size() - 1; i >= 0; i--) {
                        if ((*local_pool)[i].invalid_ptr == &(curlocal_nodes.back().node_invalid)) {
                            local_pool->erase(local_pool->begin() + i);
                            break;
                        }
                    }
                    curlocal_nodes.pop_back();
                }
            }
        }
    }
    lp_lock[thread_id].lck.unlock();

    return;
}


bool solver::grab_restartnode() {

    bool selected = false;

    if (thread_load[thread_id].next_action == abandon_action::EXPLORE_UNKNOWN) {
        while (!selected && !GPQ.Unknown.empty()) {
            if (GPQ.Unknown.back().load_info >= best_cost) {
                GPQ.Unknown.pop_back();
            }
            if (!GPQ.Unknown.empty()) {
                Generate_SolverState(GPQ.Unknown.back());
                GPQ.Unknown.pop_back();
                selected = true;
                //cout << " Thread " << thread_id << " Pick up unknown workload!. Current Initial Global Pool Size is " << GPQ.Unknown.size() << endl;
                break;
            }
            else {
                cout << "******************" << endl;
                cout << "Unknown pool is emptied" << endl;
                cout << "******************" << endl;
            }
        }
    }
    else if (thread_load[thread_id].next_action == abandon_action::EXPLOITATION_PROMISE) {
        while (!selected) {
            while (!GPQ.Promise.space.back().empty()) {
                if (GPQ.Promise.space.back().back().load_info >= best_cost || GPQ.Promise.space.back().back().deprecated) {
                    //cout << " Thread " << thread_id << " Pruned shared work!. Current Shared Pool Size is " << GPQ.Promise.space.back().size() << endl;
                    if (!GPQ.Promise.space.back().back().deprecated) {
                        GPQ.Promise.space.back().back().partial_active_path.incre_children_cnt(Allocator);
                        if (GPQ.Promise.space.back().back().root_his_node != NULL) {
                            GPQ.Promise.space.back().back().root_his_node->explored = true;
                        }
                    }
                    GPQ.Promise.space.back().pop_back();
                }
                else {
                    Generate_SolverState(GPQ.Promise.space.back().back());
                    GPQ.Promise.space.back().pop_back();
                    selected = true;
                    //cout << " Thread " << thread_id << " Pick up shared work!. Current Shared Pool Size is " << GPQ.Promise.space.back().size() << endl;
                    break;
                }
            }
            if (GPQ.Promise.space.back().empty()) {
                cout << "******************" << endl;
                cout << "Sharing pool " << GPQ.Promise.space.size() << " is emptied" << endl;
                cout << "******************" << endl;
                GPQ.Promise.space.pop_back();
                global_concentrate_lv--;
            }

            if (!selected && GPQ.Promise.is_empty()) return false;
        }
    }

    return true;
}

bool solver::Grab_from_GPQ(bool reserve) {
    int assigned = false;
    int terminate = true;
    
    GPQ_lock.lock();
    while (!assigned && !GPQ.Unknown.empty()) {
        if (GPQ.Unknown.back().load_info >= best_cost) {
            assigned = false;
        }
        else {
            Generate_SolverState(GPQ.Unknown.back());
            terminate = false;
            assigned = true;
        }
        GPQ.Unknown.pop_back();
    }

    while (!assigned && !GPQ.Abandoned.empty()) {
        if (GPQ.Abandoned.back().load_info >= best_cost || GPQ.Abandoned.back().deprecated) {
            if (!GPQ.Abandoned.back().deprecated) {
                GPQ.Abandoned.back().partial_active_path.incre_children_cnt(Allocator);
                if (GPQ.Abandoned.back().root_his_node != NULL) {
                    GPQ.Abandoned.back().root_his_node->explored = true;
                }
            }
            assigned = false;
        }
        else {
            Generate_SolverState(GPQ.Abandoned.back());
            terminate = false;
            assigned = true;
        }
        GPQ.Abandoned.pop_back();
    }
    GPQ_lock.unlock();

    
    return terminate;
}

void lkh() {
    while(!BB_Complete) {
        LKH(&f_name[0]);
        if (initial_LKHRun) {
            initial_LKHRun = false;
        }
        BB_SolFound = false;
    }
    return;
}

void solver::solve_parallel(int thread_num, int pool_size) {
    start_time_limit = std::chrono::system_clock::now();
    vector<solver> solvers(thread_num);
    deque<sop_state>* solver_container;
    vector<thread> Thread_manager(thread_num);
    vector<int> ready_thread;
    //Initially fill the GPQ with solvers:
    vector<node> ready_list;

    solver_container = new deque<sop_state>();
    problem_state.cur_solution.push_back(0);
    problem_state.taken_arr[0] = 1;
    for (int vertex : dependent_graph[0]) problem_state.depCnt[vertex]--;

    for (int i = node_count-1; i >= 0; i--) {
        if (!problem_state.depCnt[i] && !problem_state.taken_arr[i]) {
            //Push vertices with 0 depCnt into the ready list
            ready_list.push_back(node(i,-1));
        }
    }

    //Initial filling of the GPQ

    sop_state initial_state = problem_state;

    for (auto node : ready_list) {
        int taken_node = node.n;
        int cur_node = initial_state.cur_solution.back();
        initial_state.taken_arr[taken_node] = 1;
        for (int vertex : dependent_graph[taken_node]) initial_state.depCnt[vertex]--;
        initial_state.cur_cost += cost_graph[cur_node][taken_node].weight;
        initial_state.cur_solution.push_back(taken_node);
        initial_state.originate = taken_node;
        initial_state.hungarian_solver.fix_row(cur_node, taken_node);
        initial_state.hungarian_solver.fix_column(taken_node, cur_node);
        initial_state.hungarian_solver.solve_dynamic();
        initial_state.load_info = initial_state.hungarian_solver.get_matching_cost()/2;
        //number_of_lb++;
        if (pre_density != 0) solver_container->push_back(initial_state);
        else if (node_count > thread_total) {
            GPQ.Unknown.push_back(instrct_node(initial_state.cur_solution,initial_state.originate,initial_state.load_info,-1,
                                               Active_Path(initial_state.cur_solution.size()),NULL));
        }

        initial_state.taken_arr[taken_node] = 0;
        for (int vertex : dependent_graph[taken_node]) initial_state.depCnt[vertex]++;
        initial_state.cur_cost -= cost_graph[cur_node][taken_node].weight;
        initial_state.cur_solution.pop_back();
        initial_state.hungarian_solver.undue_row(cur_node, taken_node);
        initial_state.hungarian_solver.undue_column(taken_node, cur_node);
    }

    for (int i = thread_num - 1; i >= 0; i--) ready_thread.push_back(i);

    //While GPQ is not empty do split operation or assign threads with new node.
    while (GPQ.Unknown.empty() && (solver_container->size() < (size_t)pool_size || Split_level_check(solver_container))) {
        sort(solver_container->begin(),solver_container->end(),Split_sort);
        auto target = solver_container->front();

        if (target.cur_solution.size() == (unsigned) (node_count - 1)) break;
        solver_container->pop_front();
        ready_list.clear();
        for (int i = node_count-1; i >= 0; i--) {
            if (!target.depCnt[i] && !target.taken_arr[i]) ready_list.push_back(node(i,-1));
        }

        if (!ready_list.empty()) {
            for (auto node : ready_list) {
                int vertex = node.n;
                //Push split node back into GPQ
                int taken_node = vertex;
                int cur_node = target.cur_solution.back();
                target.taken_arr[taken_node] = 1;
                for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]--;
                target.cur_cost += cost_graph[cur_node][taken_node].weight;
                target.cur_solution.push_back(taken_node);
                if (target.cur_solution.size() == (size_t)node_count && target.cur_cost < best_cost) {
                    best_solution = target.cur_solution;
                    best_cost = target.cur_cost;
                }
                target.hungarian_solver.fix_row(cur_node, taken_node);
                target.hungarian_solver.fix_column(taken_node, cur_node);
                target.hungarian_solver.solve_dynamic();
                target.load_info = target.hungarian_solver.get_matching_cost()/2;
                solver_container->push_back(target);

                target.taken_arr[taken_node] = 0;
                for (int vertex : dependent_graph[taken_node]) target.depCnt[vertex]++;
                target.cur_cost -= cost_graph[cur_node][taken_node].weight;
                target.cur_solution.pop_back();
                target.hungarian_solver.undue_row(cur_node, taken_node);
                target.hungarian_solver.undue_column(taken_node, cur_node);
            }
        }
        else solver_container->push_back(target);
    }
    if (GPQ.Unknown.empty()) {
        for (auto problem : *solver_container) {
            GPQ.Unknown.push_back(instrct_node(problem.cur_solution,problem.originate,problem.load_info,-1,
                                               Active_Path(problem.cur_solution.size()),NULL));
        }
    }
    sort(GPQ.Unknown.begin(),GPQ.Unknown.end(),GPQ_sort);
    delete solver_container;

    int thread_cnt = 0;
    //cout << "GPQ initial depth is " << GPQ.back().cur_solution.size() << endl;
    cout << "Initial GPQ size is " << GPQ.Unknown.size() << endl;
    //calculate_standard_deviation();
    
    //Assign GPQ subproblem into solver;
    vector<bool> origin_taken_arr = vector<bool>(node_count,false);
    while (thread_cnt < thread_num) {
        fill(origin_taken_arr.begin(),origin_taken_arr.end(),false);
        for (int k = GPQ.Unknown.size() - 1; k >= 0; k--) {
            if (thread_cnt >= thread_num) break;
            unsigned origin = GPQ.Unknown[k].originate;
            if (!origin_taken_arr[origin]) {
                instrct_node sequence_node = GPQ.Unknown[k];
                GPQ.Unknown.erase(GPQ.Unknown.begin()+k);
                solvers[thread_cnt].problem_state = default_state;
                int taken_node = -1;
                int cur_node = sequence_node.sequence.front();
                int size = sequence_node.sequence.size();
                int counter = 0;
                solvers[thread_cnt].problem_state.cur_solution.push_back(0);
                solvers[thread_cnt].problem_state.taken_arr[cur_node] = 1;
                for (int vertex : dependent_graph[cur_node]) solvers[thread_cnt].problem_state.depCnt[vertex]--;
                counter++;

                while (counter < size) {
                    taken_node = sequence_node.sequence[counter];
                    solvers[thread_cnt].problem_state.taken_arr[taken_node] = 1;
                    for (int vertex : dependent_graph[taken_node]) solvers[thread_cnt].problem_state.depCnt[vertex]--;
                    solvers[thread_cnt].problem_state.cur_cost += cost_graph[cur_node][taken_node].weight;
                    solvers[thread_cnt].problem_state.cur_solution.push_back(taken_node);
                    solvers[thread_cnt].problem_state.hungarian_solver.fix_row(cur_node, taken_node);
                    solvers[thread_cnt].problem_state.hungarian_solver.fix_column(taken_node, cur_node);
                    solvers[thread_cnt].problem_state.initial_depth = solvers[thread_cnt].problem_state.cur_solution.size();
                    cur_node = taken_node;
                    counter++;
                }

                solvers[thread_cnt].cur_active_tree = Active_Path(solvers[thread_cnt].problem_state.cur_solution.size());
                solvers[thread_cnt].cur_active_tree.set_threadID(thread_cnt, thread_total);
                solvers[thread_cnt].problem_state.load_info = sequence_node.load_info;
                solvers[thread_cnt].lb_curlv = sequence_node.load_info;
                solvers[thread_cnt].problem_state.originate = sequence_node.originate;
                solvers[thread_cnt].node_count = node_count;
                solvers[thread_cnt].problem_state.initial_depth = solvers[thread_cnt].problem_state.cur_solution.size();
                solvers[thread_cnt].thread_id = thread_cnt;
                solvers[thread_cnt].local_pool = new deque<instrct_node>();
                thread_load[thread_cnt].cur_sol = &solvers[thread_cnt].problem_state.cur_solution;
                glp[thread_cnt].local_pool = solvers[thread_cnt].local_pool;
                origin_taken_arr[origin] = true;

                boost::dynamic_bitset<> bit_vector(node_count, false);
                for (auto node : solvers[thread_cnt].problem_state.cur_solution) {
                    bit_vector[node] = true;
                }
                int last_element = solvers[thread_cnt].problem_state.cur_solution.back();
                solvers[thread_cnt].problem_state.key = make_pair(bit_vector,last_element);

                thread_cnt++;
            }
        }
    }

    time_point = chrono::high_resolution_clock::now();

    if (enable_lkh) LKH_thread = thread(lkh);

    for (int i = 0; i < thread_num; i++) {
        Thread_manager[i] = thread(&solver::enumerate,move(solvers[i]));
        active_thread++;
    }
    
    
    for (int i = 0; i < thread_num; i++) {
        if (Thread_manager[i].joinable()) {
            Thread_manager[i].join();
            delete solvers[i].local_pool;
            ready_thread.push_back(i);
        }
    }
    active_thread = 0;

    BB_Complete = true;

    if (time_out) cout << "instance timed out " << endl;

    if (time_out || (GPQ.Unknown.empty() && GPQ.Abandoned.empty())) {
        delete thread_load;
        return;
    }

    return;
}

size_t solver::transitive_closure(vector<vector<int>>& isucc_graph) {
    size_t node_num = 0;
    deque<int> dfs_queue;
    boost::dynamic_bitset<> taken_arr = boost::dynamic_bitset<>(node_count);
    boost::dynamic_bitset<> finished_arr = boost::dynamic_bitset<>(node_count);
    vector<boost::dynamic_bitset<>> transitive_closure_map = vector<boost::dynamic_bitset<>>(node_count);
    for (int i = 0; i < node_count; i++) transitive_closure_map[i] = boost::dynamic_bitset<>(node_count);
    int selected_node = 0;

    for (int i = 0; i < node_count; i++) {
        //cout << "current node is " << i << endl;
        transitive_closure_map[i][i] = true;
        dfs_queue.push_back(i);
        while (!dfs_queue.empty()) {
            selected_node = dfs_queue.back();
            dfs_queue.pop_back();
            for (int node : isucc_graph[selected_node]) {
                if (!taken_arr[node]) {
                    taken_arr[node] = true;
                    dfs_queue.push_back(node);
                }
                if (!transitive_closure_map[i][node]) transitive_closure_map[i][node] = true;
                if (!transitive_closure_map[node][i]) transitive_closure_map[node][i] = true;
            }
        }
        taken_arr.reset();
    }

    for (int i = 0; i < node_count; i++) {
        for (int k = 0; k < node_count; k++) {
            if (transitive_closure_map[i][k]) node_num++;
        }
    }
    //print_transitive_map();
    return node_num;
}

void solver::solve(string filename,int thread_num) {
    if (thread_num == -1 || thread_num == 0) {
        cerr << "Incorrect Thread Number Input" << endl;
    }
    if (enable_lkh) thread_total = thread_num - 1;
    else thread_total = thread_num;
    f_name = filename;
    retrieve_input(filename);
    //Remove redundant edges in the cost graph
    transitive_redundantcy();
    /////////////////////////////////////////
    //Calculate Precedance Constraint Density
    int total_edges = 0;
    for (unsigned i = 1; i < dependent_graph.size() - 1; i++) {
        total_edges += dependent_graph[i].size() - 1;
    }
    pre_density = float(transitive_closure(dependent_graph)) / float((node_count)*(node_count-1));
    
    cout << "Precedance Density = " << pre_density << endl;
    /////////////////////////////////////////

    Local_PoolConfig(float((node_count - 1)*(node_count-2))/float(total_edges));

    cout << "Nearest Neighbor Heuristic Initialized" << endl;
    vector<int> temp_solution;
    temp_solution.push_back(0);
    best_solution = nearest_neightbor(&temp_solution);
    
    cout << "Initial Best Solution Is " << best_cost << endl;

    bestBB_tour = new int[node_count];
    for (int i = 0; i < node_count; i++) bestBB_tour[i] = best_solution[i] + 1;
    
    promise_Tlimit = (thread_total * exploitation_per) / tgroup_ratio;
    cout << "Maximum exploitation group during thread restart is set to " << promise_Tlimit << endl;

    max_edge_weight = get_maxedgeweight();
    problem_state.hungarian_solver = Hungarian(node_count, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    problem_state.hungarian_solver.start()/2;
    problem_state.depCnt = vector<int>(node_count,0);
    problem_state.taken_arr = vector<int>(node_count,0);

    for (int i = 0; i < node_count; i++) {
        for (unsigned k = 0; k < dependent_graph[i].size(); k++) {
            problem_state.depCnt[dependent_graph[i][k]]++;
        }
    }
    default_state = problem_state;
    cout << "Instance size is " << node_count - 2 << endl;

    thread_load = new load_stats [thread_total];

    history_table.set_up_mem(thread_total,node_count);
    
    abandon_wlk_array = vector<bool_64>(thread_total);
    num_resume = vector<int_64>(thread_total);
    num_stop = vector<int_64>(thread_total);
    glp = vector<lptr_64>(thread_total);
    lp_lock = vector<mutex_64>(thread_total);
    busy_arr = vector<atomwrapper<bool>>(thread_total);

    steal_wait = vector<double>(thread_total,0);
    lp_time = vector<double>(thread_total,0);
    proc_time = vector<vector<double>>(thread_total);
    for (int i = 0; i < thread_total; i++) proc_time[i] = vector<double>(3,0);
    for (int i = 0; i < thread_total; i++) initial_hungstate.push_back(default_state.hungarian_solver);
    steal_cnt = vector<int>(thread_total,0);
    enumerated_nodes = vector<unsigned_long_64>(thread_total);
    
    auto start_time = chrono::high_resolution_clock::now();
    solve_parallel(thread_total,global_pool_size);
    auto end_time = chrono::high_resolution_clock::now();

    if (enable_lkh) if (LKH_thread.joinable()) LKH_thread.join();

    auto total_time = chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    cout << "------------------------" << thread_total << " thread" << "------------------------------" << endl;
    cout << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;
    
    
    /*
    for (int i = 0; i < thread_total; i++) {
        cout << "--------------------------------------------------------" << endl;
        cout << "Thread " << i << ": Steal Count = " << steal_cnt[i] << endl;
        //cout << "                    Wait Time = " << steal_wait[i] / (float)(1000000) << " s" << endl;
        //cout << "                    Local Pool Time = " << lp_time[i] / (float)(1000000) << " s" << endl;
        //cout << "                    Load Process Time 0 = " << proc_time[i][0] / (float)(1000000) << " s" << endl;
        //cout << "                    Load Process Time 1 = " << proc_time[i][1] / (float)(1000000) << " s" << endl;
        //cout << "                    Load Process Time 2 = " << proc_time[i][2] / (float)(1000000) << " s" << endl;
        //cout << "--------------------------------------------------------" << endl;
    }
    */
    
    /*
    unsigned long total_enodes = 0;
    for (int i = 0; i < thread_total; i++) total_enodes += enumerated_nodes[i].val;

    for (int i = 0; i < thread_total; i++) {
        cout << "thread " << i << " enumerated nodes = " << enumerated_nodes[i].val << ", " << double(enumerated_nodes[i].val)/double(total_enodes) << "%" << endl;
    }
    */
    

    //history_table.print_curmem();
    //cout << "Total enumerated nodes are " << total_enodes << endl;
    /*
    if (best_cost != 71749) {    
        for (int i = 0; i < thread_total; i++) {
            cout << "Thread id == " << i << endl;
            for (int k = 0; k < tree[i].size(); k++) {
                if (tree[i][k].HistoryPruned) cout << "HX: cc = " 
                                                   << tree[i][k].cur_cost 
                                                   << ", pc = " << tree[i][k].prefix_cost 
                                                   << ", hlb = " << tree[i][k].hlb
                                                   << ", best_cost = " << tree[i][k].best_cost
                                                   << ", athread = " << tree[i][k].active_thread << ": ";
                else if (tree[i][k].LBPruned) cout << "LBX plb = " << tree[i][k].plb << ", best_cost = " << tree[i][k].best_cost << ": ";
                else cout << "pc = " << tree[i][k].prefix_cost << ", lb = " << tree[i][k].plb << " :";
                for (auto node : tree[i][k].solution) cout << node << ",";
                cout << endl;
            }
            cout << endl;
        }
    }
    */
        
    return;
}


void solver::Local_PoolConfig(float precedence_density) {
    int min_lpn = 100;
    if (pre_density < 0.5) {
        if (node_count < 500) local_pool_size = min_lpn * 1.5;
        else local_pool_size = min_lpn * 2.5;
    }
    else {
        if (node_count < 100) local_pool_size = min_lpn * 1.5;
        else if (node_count < 500) local_pool_size = min_lpn * 3;
        else local_pool_size = min_lpn * 4;
    }
    
    if (!enable_workstealing) {
        local_pool_size = 0;
        cout << "Work stealing disabled" << endl;
    }
    else cout << "Work stealing enabled" << endl;
    
    if (!enable_threadstop) cout << "Thread stopping disabled" << endl;
    else cout << "Thread stopping enabled" << endl;

    cout << "LPQ size = " << local_pool_size << endl;

    return;
}

void solver::retrieve_input(string filename) {
    ifstream inFile;
    string line;
    string fd_name;

    for (unsigned i = 0; i < filename.size(); i++) {
        if (i == filename.size() - 4) break;
        fd_name += filename[i];
    }
    fd_name += "_BB.sop";

    inFile.open(fd_name);

    if (inFile.fail()) {
        cerr << "Error: input file " << fd_name << " -> " << strerror(errno) << endl;
        exit(-1);
    }
    else cout << "Input file is " << fd_name << endl;

    // Read input files and store it inside an array.
    vector<vector<int>> file_matrix;
    vector<int> matrix_temp;
    while (getline(inFile,line)) {  
        stringstream sstream;
        sstream << line;
        string weight;
        int weight_num;
        while (sstream >> weight) {
            stringstream(weight) >> weight_num;
            matrix_temp.push_back(weight_num);
        }
        file_matrix.push_back(matrix_temp);
        matrix_temp.clear();
    }

    unsigned size = file_matrix.size();
    
    cost_graph = vector<vector<edge>>(size);
    hung_graph = vector<vector<edge>>(size);
    dependent_graph = vector<vector<int>>(size);
    
    for (int i = 0; i < (int)file_matrix.size(); i++) {
        int j = 0;
        for (auto edge_weight: file_matrix[i]) {
            if (edge_weight < 0) {
                cost_graph[i].push_back(edge(i,j,file_matrix[j][i]));
                hung_graph[i].push_back(edge(i,j,-1));
                dependent_graph[j].push_back(i);
            }
            else { 
                cost_graph[i].push_back(edge(i,j,edge_weight));
                hung_graph[i].push_back(edge(i,j,edge_weight));
            }
            j++;
        }
    }
    node_count = cost_graph.size();
    //Trim redundant edges
    return;
}

void solver::transitive_redundantcy() {
    in_degree = std::vector<vector<edge>>(node_count);
    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            int c = dependent_graph[i][k];
            in_degree[c].push_back(edge(i,c,hung_graph[i][c].weight));
        }
    }

    for(int i = 0; i < node_count; ++i) {
        vector<edge> preceding_nodes;
        for (int k = 0; k < (int)dependent_graph[i].size(); k++) preceding_nodes.push_back(edge(i,dependent_graph[i][k],-1));
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hung_graph[dependence_edge.dest][i].weight = -1;
                    hung_graph[i][dependence_edge.dest].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }

                for(int dest : dependent_graph[dependence_edge.dest]){
                    if(expanded_nodes.find(dest) == expanded_nodes.end()){
                        st.push_back(edge(dependence_edge.dest,dest,-1));
                        expanded_nodes.insert(dest);
                    }
                }
            }
        } 
    }

    for(int i = 0; i < node_count; ++i) {
        const vector<edge> preceding_nodes = in_degree[i];
        unordered_set<int> expanded_nodes;
        for(int j = 0; j < (int)preceding_nodes.size(); ++j){
            vector<edge> st;
            st.push_back(preceding_nodes[j]);
            while(!st.empty()){
                edge dependence_edge = st.back();
                st.pop_back();
                if(dependence_edge.src != i){
                    hung_graph[i][dependence_edge.dest].weight = -1;
                    hung_graph[dependence_edge.dest][i].weight = -1;
                    expanded_nodes.insert(dependence_edge.dest);
                }
                for(const edge& e : in_degree[dependence_edge.dest]){
                    if(expanded_nodes.find(e.dest) == expanded_nodes.end()){
                        st.push_back(e);
                        expanded_nodes.insert(e.dest);
                    }
                }
            }
        } 
    }

    return;
}


bool compare (const edge src, const edge target) {
    int src_weight = src.weight;
    int dest_weight = target.weight;
    return (src_weight < dest_weight);
}

//Sort cost-graph weight by ascending order
void solver::sort_weight(vector<vector<edge>>& graph) {
    int size = graph.size();
    for (int i = 0; i < size; i++) {
        stable_sort(graph[i].begin(),graph[i].end(),compare);
    }
    return;
}


vector<int> solver::nearest_neightbor(vector<int>* partial_solution) {
    vector<int> solution;
    int current_node;
    int solution_cost = 0;
    bool visit_arr[node_count];
    int depCnt_arr[node_count];
    vector<vector<edge>> sorted_costgraph = cost_graph; 
    sort_weight(sorted_costgraph);
    //sort input based on weight for NN heurestic
    memset(depCnt_arr,0,node_count*sizeof(int));

    for (int i = 0; i < node_count; i++) {
        visit_arr[i] = false;
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    current_node = 0;
    for (auto node : *partial_solution) {
        visit_arr[node] = true;
        for (long unsigned int i = 0; i < dependent_graph[node].size(); i++) {
            depCnt_arr[dependent_graph[node][i]]--;
        }
        solution.push_back(node);
        solution_cost += cost_graph[current_node][node].weight;
        current_node = node;
    }

    int num = solution.size();
    
    while (num < node_count) {
        bool taken = false;
        for (auto node: sorted_costgraph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest]) {
                current_node = node.dest;
                solution_cost += node.weight;
                solution.push_back(current_node);
                num++;
                visit_arr[node.dest] = true;
                taken = true;
                break;
            }
        }
        if (!taken) {
            cout << "current node is " << current_node << endl;
            for (auto node: sorted_costgraph[current_node]) {
                cout << node.dest << "," << visit_arr[node.dest] << "," << depCnt_arr[node.dest] << endl;
            }
            exit(1);
        }
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
    }

    best_cost = solution_cost;
    return solution;
}


int solver::get_maxedgeweight() {
    int max = 0;
    for (int i = 0; i < node_count; i++) {
        for (auto edge_weight : cost_graph[i]) {
            int weight = edge_weight.weight;
            if (weight > max) max = weight;
        }
    }
    return max;
}

vector<vector<int>> solver::get_cost_matrix(int max_edge_weight) {
    vector<vector<int>> matrix(node_count);
    for(int i = 0; i < node_count; ++i){
		matrix[i] = vector<int>(node_count, max_edge_weight*2);
	}

    for (vector<edge> edge_list : hung_graph) {
        for (auto edge : edge_list) {
            int i = edge.src;
            int k = edge.dest;
            int weight = edge.weight;
            if (weight != -1 && i != k) {
                matrix[i][k] = weight * 2;
            }
        }
    }

    return matrix;
}