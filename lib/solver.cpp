#include "solver.hpp"
#include "hash.hpp"
#include <iostream>
#include <fstream>
#include <stack>
#include <sstream>
#include <limits>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <functional>
#include <vector>
#include <bits/stdc++.h>
#include <unordered_map>
#include <setjmp.h>

using namespace std;
#define TABLE_SIZE 541065431
#define BLOCK_SIZE 81920
#define ENABLE_PROGRESS_ESTIMATION 1            //tracks what percentage of the working tree has been traversed or pruned so far
#define PROGRESS_ESTIMATION_USE_REMAINDER 0     //if progress estimation is enabled, will use the more accurate but marginally more expensive algorithm

Hash_Map history_table(TABLE_SIZE);
bool optimal_found = false;
int graph_index = 0;
int recently_added = 0;
long eclipsed_time = 0;

HistoryNode visted;
bool full_solution = false;
float total_time_static = 0;
float ready_list_time = 0;
int lb = 0;
int suffix_cost = 0;
int previous_snode = 0;

std::chrono::time_point<std::chrono::system_clock> start_time_limit;
jmp_buf buf;
vector<int> raedy_list;
vector<int> suffix;
vector<bool> picked_list;

float wait_time = 0;
unsigned long enumerated_nodes = 0;

#if ENABLE_PROGRESS_ESTIMATION
    unsigned long long current_node_value = ULLONG_MAX; //the portion out of ULLONG_MAX of the working tree that is under each child
    unsigned long long leaf_percent = 0; //total portion of the tree, out of ULLONG_MAX, that was fully processed
    unsigned long long trimmed_backtracking = 0; //accumulated estimated portion of entire tree trimmed, not including leaves, out of ULLONG_MAX
    unsigned long long trimmed_history = 0; //accumulated estimated portion of entire tree trimmed, not including leaves, out of ULLONG_MAX
    unsigned long long trimmed_hungarian = 0; //accumulated estimated portion of entire tree trimmed, not including leaves, out of ULLONG_MAX
    #if PROGRESS_ESTIMATION_USE_REMAINDER
        bool backtracking_complete = false, history_complete = false, hungarian_complete = false;
    #endif
#endif

int* depCnt;
int* taken_arr;

void solver::push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,int i) {
    if (history_table.get_cur_size() < 0.8 * history_table.get_max_size()) history_table.insert(key,cur_cost,lower_bound,best_cost);
    else if (i < int(0.5 * node_count) && history_table.get_cur_size() < history_table.get_max_size()) history_table.insert(key,cur_cost,lower_bound,best_cost);
    return;
}

int solver::dynamic_hungarian(int src, int dest) {
    hungarian_solver.fix_row(src, dest);
	hungarian_solver.fix_column(dest, src);
	hungarian_solver.solve_dynamic();
    return hungarian_solver.get_matching_cost()/2;
}


bool solver::HistoryUtilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,int cost) {
    size_t val = hash<boost::dynamic_bitset<>>{}(key.first);
    int bucket_location = (val + key.second) % TABLE_SIZE;
    HistoryNode* history_node = history_table.retrieve(key,bucket_location);

    if (history_node == NULL) return true;

    *found = true;
    HistoryContent Loaded_Val = history_node->Entry;
    int history_prefix = Loaded_Val.prefix_cost;
    int history_lb = Loaded_Val.lower_bound;
    *lowerbound = history_lb;

    if (cost >= history_prefix) return false;
    
    int imp = history_prefix - cost;
    
    if (imp <= history_lb - best_cost) return false;

    history_node->Entry = {cost,history_lb - imp};
    *lowerbound = history_lb - imp;

    return true;
}

int solver::dynamic_edb() {
    int picked_node = cur_solution.back();
    picked_list[picked_node] = true;
    recently_added = picked_node;
    return mmcp_lb();
}

bool bound_sort(const node& src,const node& dest) {
    if (src.lb == dest.lb) return src.nc > dest.nc;
    return src.lb > dest.lb;
}

bool nearest_sort(const node& src,const node& dest) {
    if (src.nc == dest.nc) return src.lb > dest.lb;
    return src.nc > dest.nc;
}

void solver::enumerate(int depth) {
    vector<node> ready_list;

    for (int i = node_count-1; i >= 0; i--) {
        if (!depCnt[i] && !taken_arr[i]) {
            ready_list.push_back(node(i,-1,-1));
        }
    }

    #if ENABLE_PROGRESS_ESTIMATION
        unsigned long long parent_node_value = current_node_value; //the portion out of ULLONG_MAX of the working tree that is under this node, which will be split between its children
        current_node_value = parent_node_value / ready_list.size(); //the parent's value is assumed to be split evenly between all children, even though this is not strictly correct
        #if PROGRESS_ESTIMATION_USE_REMAINDER
            int remaining_share = parent_node_value % ready_list.size(); //remainder never makes enough of a difference to even be visible
        #endif
        //unsigned long long local_progress_accumulator = 0; //simply for testing if the progress estimation is working correctly
    #endif

    int last_element = cur_solution.back();
    int taken_node = 0;
    int u = 0;
    int v = 0;

    deque<node> enumeration_list;

    enumerated_nodes += ready_list.size();
    #if ENABLE_PROGRESS_ESTIMATION
        #if PROGRESS_ESTIMATION_USE_REMAINDER
            int remainder_not_pruned = 0; //number of children from the first X % Y children that is added to the enumeration_list
        #endif
    #endif

    if (!ready_list.empty()) { //then do pruning
        for (int i = 0; i < (int)ready_list.size(); i++) {
            node dest = ready_list[i];
            int src = cur_solution.back();
            cur_solution.push_back(dest.n);
            cur_cost += cost_graph[src][dest.n].weight;
            int temp_lb = -1;
            bool taken = false;

            #if ENABLE_PROGRESS_ESTIMATION
                #if PROGRESS_ESTIMATION_USE_REMAINDER
                    if (remaining_share > 0) //if this node's share of the solution space cannot be evenly split between children
                    {
                        if (i == 0) //the first X % Y children have their value increased by 1
                            current_node_value++;
                        else if (i == remaining_share) //the remaining children have a slightly lower value arbitrarily, guaranteeing that no share is lost by rounding
                            current_node_value--;
                    }
                #endif
            #endif
            
            if (cur_cost >= best_cost) {
                #if ENABLE_PROGRESS_ESTIMATION
                    #if PROGRESS_ESTIMATION_USE_REMAINDER
                        unsigned long long tmp = trimmed_backtracking + current_node_value;
                        if (tmp < trimmed_backtracking)
                            backtracking_complete = true;
                        trimmed_backtracking = tmp;
                    #else
                        trimmed_backtracking += current_node_value;
                    #endif
                    //local_progress_accumulator += current_node_value;
                #endif
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                continue;
            }
            if (cur_solution.size() == (unsigned)node_count) {
                if (cur_cost < best_cost) {
                    best_solution = cur_solution;
                    best_cost = cur_cost;
                    cout << "best cost = " << best_cost << " found at time = " << std::chrono::duration<double>(std::chrono::system_clock::now() - start_time_limit).count() << " seconds " << endl;
                }
                #if ENABLE_PROGRESS_ESTIMATION
                    leaf_percent += current_node_value; 
                    //local_progress_accumulator += current_node_value;
                #endif
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                continue;
            }
            
            else {
                key.first[dest.n] = true;
                key.second = dest.n;
                bool decision = HistoryUtilization(key,&temp_lb,&taken,cur_cost);
                if (!taken) { //if there is no similar entry in the history table
                    temp_lb = dynamic_hungarian(src,dest.n);
                    push_to_historytable(key,temp_lb,i);
                    hungarian_solver.undue_row(src,dest.n);
                    hungarian_solver.undue_column(dest.n,src);
                }
                else if (taken && !decision) { //if this path is dominated by another path
                    #if ENABLE_PROGRESS_ESTIMATION
                        #if PROGRESS_ESTIMATION_USE_REMAINDER
                            unsigned long long tmp = trimmed_history + current_node_value;
                            if (tmp < trimmed_history)
                                history_complete = true;
                            trimmed_history = tmp;
                        #else
                            trimmed_history += current_node_value;
                        #endif
                        //local_progress_accumulator += current_node_value;
                    #endif
                    key.first[dest.n] = false;
                    key.second = last_element;
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    continue;
                }
                if (temp_lb >= best_cost) { //pruning based on dynamic hungarian (technically overestimates because history domination can find a better lower bound)
                    #if ENABLE_PROGRESS_ESTIMATION
                        #if PROGRESS_ESTIMATION_USE_REMAINDER
                            unsigned long long tmp = trimmed_hungarian + current_node_value;
                            if (tmp < trimmed_hungarian)
                                hungarian_complete = true;
                            trimmed_hungarian = tmp;
                        #else
                            trimmed_hungarian += current_node_value;
                        #endif
                        //local_progress_accumulator += current_node_value;
                    #endif
                    cur_solution.pop_back();
                    cur_cost -= cost_graph[src][dest.n].weight;
                    key.first[dest.n] = false;
                    key.second = last_element;
                    continue;
                }

                //if you don't prune a node, add it to the enumeration_list
                cur_solution.pop_back();
                cur_cost -= cost_graph[src][dest.n].weight;
                enumeration_list.push_back(ready_list[i]);
                enumeration_list.back().nc = cost_graph[src][dest.n].weight;
                enumeration_list.back().lb = temp_lb;
                key.first[dest.n] = false;
                key.second = last_element;
                #if ENABLE_PROGRESS_ESTIMATION
                    #if PROGRESS_ESTIMATION_USE_REMAINDER
                        if (remaining_share > 0 && i < remaining_share) 
                            remainder_not_pruned++;
                    #endif
                #endif
            }
        }
        sort(enumeration_list.begin(),enumeration_list.end(),bound_sort);
    }

    // #if ENABLE_PROGRESS_ESTIMATION
    //     unsigned long long should_be = local_progress_accumulator + enumeration_list.size()*current_node_value + remainder_not_pruned;
    //     if (should_be != parent_node_value)
    //     {
    //         cout << "PROGRESS LEAK. " << should_be << " != " << parent_node_value << endl;
    //         cout << "Current path: ";
    //         for (int i : cur_solution)
    //         {
    //             cout << i << " ";
    //         }
    //         cout << endl;
    //         cout << "Parent value: " << parent_node_value << endl;
    //         cout << "Pruned: " << local_progress_accumulator << endl;
    //         cout << "Size of enumeration list: " << enumeration_list.size() << endl;
    //         cout << "Child value (floor): " << current_node_value << endl;
    //         cout << "Remainder: " << remainder_not_pruned << endl;
    //     }
    // #endif

    int child_num = 0; //the order of each child in this node's enumeration list, 0-indexed, to determine how to split inherited solution space
    while(!enumeration_list.empty()) {

        #if ENABLE_PROGRESS_ESTIMATION
            #if PROGRESS_ESTIMATION_USE_REMAINDER
                if (remainder_not_pruned > 0) //if this node's share of the solution space cannot be evenly split between children
                {
                    if (child_num == 0) //the first X % Y children have their value increased by 1
                        current_node_value++;
                    else if (child_num == remainder_not_pruned) //the remaining children have a slightly lower value arbitrarily, guaranteeing that no share is lost by rounding
                        current_node_value--;
                }
            #endif
        #endif
        
        if (enumeration_list.back().lb >= best_cost) { //you have to check again because best_cost may have changed
            #if ENABLE_PROGRESS_ESTIMATION
                #if PROGRESS_ESTIMATION_USE_REMAINDER
                    unsigned long long tmp = trimmed_hungarian + current_node_value;
                    if (tmp < trimmed_hungarian)
                        hungarian_complete = true;
                    trimmed_hungarian = tmp;
                #else
                    trimmed_hungarian += current_node_value;
                #endif
                //local_progress_accumulator += current_node_value;
            #endif
            enumeration_list.pop_back();
            continue;
        }

        taken_node = enumeration_list.back().n;
        key.first[taken_node] = true;
        key.second = taken_node;
        enumeration_list.pop_back();
        u = cur_solution.back();
        v = taken_node;
        hungarian_solver.fix_row(u, v);
        hungarian_solver.fix_column(v, u);
        cur_cost += cost_graph[u][v].weight;
        

        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]--;
        cur_solution.push_back(taken_node);
        taken_arr[taken_node] = 1;
        full_solution = false;
        suffix_cost = 0;

        enumerate(depth+1);

        for (int vertex : dependent_graph[taken_node]) depCnt[vertex]++;
        taken_arr[taken_node] = 0;
        cur_solution.pop_back();
        cur_cost -= cost_graph[u][v].weight;
        hungarian_solver.undue_row(u,v);
        hungarian_solver.undue_column(v,u);
        key.first[taken_node] = false;
        key.second = cur_solution.back();

        auto cur_time = std::chrono::system_clock::now();
        if (std::chrono::duration<double>(cur_time - start_time_limit).count() > t_limit) {
            longjmp(buf, 1);
        }
        
       // cout << "assign to history table time: " << setprecision(4) << total_time / (float)(1000000) << endl;
    }
    #if ENABLE_PROGRESS_ESTIMATION
        current_node_value = parent_node_value;
    #endif
    child_num++;
    return;
}


void solver::solve(string filename,long time_limit) {
    retrieve_input(filename);
    t_limit = time_limit;
    best_solution = nearest_neightbor();
    int max_edge_weight = get_maxedgeweight();
    hungarian_solver = Hungarian(node_count, max_edge_weight+1, get_cost_matrix(max_edge_weight+1));
    MMCP_static_lowerbound = hungarian_solver.start()/2;
    depCnt = new int[node_count];
    taken_arr = new int[node_count];
    memset(depCnt,0,node_count * sizeof(int));
    memset(taken_arr,0,node_count * sizeof(int));

    for (int i = 0; i < node_count; i++) {
        for (unsigned k = 0; k < dependent_graph[i].size(); k++) {
            depCnt[dependent_graph[i][k]]++;
        }
    }

    cur_solution.push_back(0);
    taken_arr[0] = 1;
    for (int vertex : dependent_graph[0]) depCnt[vertex]--;

    boost::dynamic_bitset<> bit_vector(node_count);
    for (auto node : cur_solution) {
        bit_vector[node] = true;
    }
    int last_element = cur_solution.back();
    key = make_pair(bit_vector,last_element);
    /*
    cout << "best solution found using NN is " << best_cost << endl;
    cout << "the NN solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    cout << "Edge-based LB is " << EGB_static_lowerbound << endl;
    cout << "MMCP-based LB is " << MMCP_static_lowerbound << endl;
    */
    cur_cost = 0;
    cout << "best initial cost is " << best_cost << endl;
    auto start_time = chrono::high_resolution_clock::now();
    start_time_limit = std::chrono::system_clock::now();
    if (setjmp(buf)) {
        cout << "Timed out!" << endl;
    }
    else {
        enumerate(1);
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<std::chrono::microseconds>( end_time - start_time ).count();

    cout << best_cost << "," << setprecision(4) << total_time / (float)(1000000) << endl;

    cout << "Total enumerated nodes are " << enumerated_nodes << endl;

    #if ENABLE_PROGRESS_ESTIMATION
        unsigned long long estimated_trimmed_percent = trimmed_backtracking + trimmed_history + trimmed_hungarian;
        unsigned long long total_progress = estimated_trimmed_percent + leaf_percent;
        #if PROGRESS_ESTIMATION_USE_REMAINDER
            bool trim_complete = estimated_trimmed_percent < trimmed_backtracking || estimated_trimmed_percent < trimmed_history || estimated_trimmed_percent < trimmed_hungarian
                                    || backtracking_complete || history_complete || hungarian_complete;
            bool final_complete = total_progress < estimated_trimmed_percent || total_progress < leaf_percent
                                    || trim_complete;
            cout << "Total Enumerated Percent = "       << ((double) leaf_percent)/ULLONG_MAX*100 << "%" << endl;
            cout << "Total Trimmed Percent = "          << (trim_complete ? 100 : ((double) estimated_trimmed_percent)/ULLONG_MAX*100) << "%" << endl;
            cout << "Trimmed in Backtracking = "        << (backtracking_complete ? 100 : ((double) trimmed_backtracking)/ULLONG_MAX*100) << "%" << endl;
            cout << "Trimmed from History Table = "     << (history_complete ? 100 : ((double) trimmed_history)/ULLONG_MAX*100) << "%" << endl;
            cout << "Trimmed by Hungarian Algorithm = " << (hungarian_complete ? 100 : ((double) trimmed_hungarian)/ULLONG_MAX*100) << "%" << endl;
            cout << "Total Progress = "                 << (final_complete ? 100 : ((double) total_progress)/ULLONG_MAX*100) << "%" << endl;
            cout << estimated_trimmed_percent + leaf_percent << " / " << ULLONG_MAX << endl;
        #else
            cout << "Total Enumerated Percent = "       << ((double) leaf_percent)/ULLONG_MAX*100 << "%" << endl;
            cout << "Total Trimmed Percent = "          << ((double) estimated_trimmed_percent)/ULLONG_MAX*100 << "%" << endl;
            cout << "Trimmed in Backtracking = "        << ((double) trimmed_backtracking)/ULLONG_MAX*100 << "%" << endl;
            cout << "Trimmed from History Table = "     << ((double) trimmed_history)/ULLONG_MAX*100 << "%" << endl;
            cout << "Trimmed by Hungarian Algorithm = " << ((double) trimmed_hungarian)/ULLONG_MAX*100 << "%" << endl;
            cout << "Total Progress = "                 << ((double) total_progress)/ULLONG_MAX*100 << "%" << endl;
            cout << estimated_trimmed_percent + leaf_percent << " / " << ULLONG_MAX << endl;
        #endif
    #endif

    //history_table.average_size();
    //cout << "Total Time Spent On HistoryUtilization Calls is " << wait_time << " s" << endl;
    //cout << "Total Number Of HistoryUtilization Calls is " << num_of_hiscall << endl;
    //cout << "Average History Utilization Call Time is " << float(wait_time / num_of_hiscall) << endl;

    /*
    cout << "the optimal solution contains ";
    for (int i = 0; i < node_count; i++) {
        if (i != node_count - 1) cout << best_solution[i] << "-->";
        else cout << best_solution[i];
    }
    cout << endl;
    */

    
    //cout << "Total number of history table insertions are " << his_insertion << endl;
    //history_table.average_size();
    //cout << "Table look up: " << number_of_history << ", LB calculation: " << number_of_lb << endl;
}

void solver::retrieve_input(string filename) {
    ifstream inFile;
    string line;
    inFile.open(filename);

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
    transitive_redundantcy();
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

    for(int i = 0; i < node_count; ++i){
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
        graph_index = i;
        stable_sort(graph[i].begin(),graph[i].end(),compare);
    }

    return;
}

vector<int> solver::nearest_neightbor() {
    vector<int> solution;
    int current_node;
    bool visit_arr[node_count];
    bool selected = false;
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

    for (int i = 0; i < node_count; i++) {
        if (depCnt_arr[i] == 0 && !selected) {
            current_node = i;
            visit_arr[current_node] = true;
            break;
        }
    }

    solution.push_back(current_node);
    
    int num = 1;
    int solution_cost = 0;
    
    while (num < node_count) {
        for (long unsigned int i = 0; i < dependent_graph[current_node].size(); i++) {
            depCnt_arr[dependent_graph[current_node][i]]--;
        }
        for (auto node: sorted_costgraph[current_node]) {
            if (!visit_arr[node.dest] && !depCnt_arr[node.dest]) {
                current_node = node.dest;
                solution_cost += node.weight;
                solution.push_back(current_node);
                num++;
                visit_arr[node.dest] = true;
                break;
            }
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

void solver::print_dep() {
    for (int i = 0; i < node_count; i++) {
        cout << "node " << i << "has children: ";
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            if (k != dependent_graph[i].size() - 1) cout << dependent_graph[i][k] << ",";
            else cout << dependent_graph[i][k];
        }
        cout << endl;
    }
}


int solver::mmcp_lb() {
    int outsum = 0;
    int insum = 0;
    int out_array[node_count];
    int in_array[node_count];
    int depCnt_arr[node_count];
    bool trim_in = true;
    bool trim_out = true;
    
    memset(depCnt_arr,0,node_count*sizeof(int));
    for (int i = 0; i < node_count; i++) {
        for (long unsigned int k = 0; k < dependent_graph[i].size(); k++) {
            depCnt_arr[dependent_graph[i][k]]++;
        }
    }

    int out_max = 0;
    int in_max = 0;

    out_array[0] = 0;
    out_array[node_count-1] = 0;
    in_array[0] = 0;
    in_array[node_count-1] = 0;


    for (int i = 1; i < node_count - 1; i++) {
        int min_out = INT_MAX;
        int min_in = INT_MAX;
        bool calculate_in = true;
        bool calculate_out = true;


        if (picked_list[i]) {
            calculate_in = false;
            trim_in = false;
            if (recently_added != i) {
                calculate_out = false;
                trim_out = false;
            }
        }


        for (int k = 1; k < node_count - 1; k++) {
            int weight = cost_graph[i][k].weight;
            if (calculate_out) {
                if (find(in_degree[i].begin(),in_degree[i].end(),edge(k,i,-1)) == in_degree[i].end()) {
                    if (k != i && weight < min_out && !picked_list[k]) min_out = weight;
                }
            }
            
            if (calculate_in) {
                if (find(in_degree[k].begin(),in_degree[k].end(),edge(i,k,-1)) == in_degree[k].end()) {
                    if (recently_added == k && k != i && weight < min_in) {
                        min_in = weight;
                    }
                    else if (!picked_list[k] && k != i && weight < min_in) {
                        min_in = weight;
                    }
                }
            }
            
        }

        if (min_out == INT_MAX) {
            trim_out = false;
            min_out = 0;
        }

        if (min_in == INT_MAX) {
            trim_in = false;
            min_out = 0;
        }

        
        if (calculate_in) in_array[i] = min_in;
        else in_array[i] = 0;

        if (calculate_out) out_array[i] = min_out;
        else out_array[i] = 0;
        
    }

    for (int i = 1; i < node_count - 1; i++) {
        if (dependent_graph[i].size() == 1 && dependent_graph[i][0] == node_count-1 && out_array[i] > out_max) {
            out_max = out_array[i];
        }
        if (depCnt_arr[i] == 1 && in_degree[i][0].src == 0 && in_array[i] > in_max) {
            in_max = in_array[i];
        }
        outsum += out_array[i];
        insum += in_array[i];
    }

    if (trim_in == false && trim_out == false) return max(outsum,insum);
    else if (trim_in == false) return max(outsum-out_max,insum);
    else if (trim_out == false) return max(outsum,insum-in_max);
    
    return max(outsum-out_max,insum-in_max);
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
