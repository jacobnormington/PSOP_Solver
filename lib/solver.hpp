#ifndef SOLVER_H
#define SOLVER_H

#include <map>
#include <string>
#include <vector>
#include <chrono>
#include <boost/dynamic_bitset.hpp>
#include "history.hpp"
#include "../lib/hungarian.h"

using namespace std;

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
        int n; //node number
        int lb; //lower bound cost of this node
        int nc; //cost to get from the parent to this node
        node(int x, int y, int z): n{x},lb{y},nc{z} {}
};

class solver {
    private:
        string enum_option;
        int node_count;
        int EGB_static_lowerbound;
        int MMCP_static_lowerbound;
        int cur_cost;
        int best_cost;
        long t_limit;

        vector<int> best_solution;
        vector<int> cur_solution;
        vector<vector<int>> dependent_graph;
        vector<vector<edge>> in_degree;
        vector<vector<edge>> hung_graph;
        vector<vector<edge>> cost_graph;
        pair<boost::dynamic_bitset<>,int> key;
        Hungarian hungarian_solver;
    public:
        vector<vector<int>> get_cost_matrix(int max_edge_weight);
        vector<int> nearest_neightbor();
        void enumerate(int i);
        bool LB_Check(int src, int dest);
        bool HistoryUtilization(pair<boost::dynamic_bitset<>,int>& key,int* lowerbound,bool* found,int cost);
        int get_maxedgeweight();
        int dynamic_edb();
        int dynamic_hungarian(int src, int dest);
        int mmcp_lb();
        int History_LB();
        void process_solution();
        void push_to_historytable(pair<boost::dynamic_bitset<>,int>& key,int lower_bound,int i);
        void assign_historytable(int prefix_cost,int lower_bound,int i);
        void solve(string filename,long time_limit);
        void retrieve_input(string filename);
        void transitive_redundantcy();
        void sort_weight(vector<vector<edge>>& graph);
        void print_dep();
        //double get_estimated_trimmed_percent(int node_count, int depth);
};

#endif