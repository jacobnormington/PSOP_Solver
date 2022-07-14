#ifndef HASH_H
#define HASH_H

#include <vector>
#include <unordered_map>
#include <list>
#include <boost/dynamic_bitset.hpp>
#include "history.hpp"

#define BUCKET_BLK_SIZE 81920
#define HIS_BLK_SIZE 81920
#define MEMORY_RESTRIC 0.85
typedef list<pair<pair<boost::dynamic_bitset<>,int>,HistoryNode*>> Bucket;

class Mem_Allocator {
    private:
        Bucket* Bucket_blk = NULL;
        unsigned counter;
        HistoryNode* history_block = NULL;
        unsigned node_counter = 0;
        int padding[8];
    public:
        Mem_Allocator();
        Bucket* Get_bucket();
        HistoryNode* retrieve_his_node();
};


class Hash_Map {
    private:
        size_t cur_size;
        size_t size = 0;
        size_t max_size = 0;
        size_t node_size = 0;
        unsigned maximum_usage_cnt = 0;
        vector<Bucket*> History_table;
        Mem_Allocator Mem_Manager;
    public:
        Hash_Map(int size);
        bool isEmpty();
        size_t get_max_size();
        size_t get_cur_size();
        uint32_t hash_func(pair<boost::dynamic_bitset<>,int> item);
        bool insert(pair<boost::dynamic_bitset<>,int>& item,int prefix_cost,int lower_bound, int best_cost);
        void increase_size(size_t size_incre);
        void set_node_t(int node_size);
        void set_up_mem(int thread_num);
        void average_size();
        void analyze_usage();
        HistoryNode* retrieve(pair<boost::dynamic_bitset<>,int>& item,int key);
};

#endif
