#ifndef HASH_H
#define HASH_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <mutex>
#include <cmath>
#include <climits>
#include <cerrno>
#include <atomic>
#include <sys/sysinfo.h>
#include <pthread.h>
#include <boost/smart_ptr/detail/spinlock.hpp>
#include <boost/dynamic_bitset.hpp>
#include "history.hpp"

#define BUCKET_BLK_SIZE 81920
#define HIS_BLK_SIZE 81920
#define FREED_SIZE 100000000
#define COVER_AREA 10
#define DEPLETION_TIME_FRAME 10

typedef pair<pair<boost::dynamic_bitset<>,int>,HistoryNode*> Entry;
typedef list<Entry> Bucket;
static mutex Gretrive_Lock;
static atomic<int> GlobalRecycle_Size(0);
static std::chrono::time_point<std::chrono::system_clock> start_clean_limit;
static std::chrono::time_point<std::chrono::system_clock> start_detect_limit;

struct Recycled_Blk {
    HistoryNode** Node_Blk = NULL;
    int total_counter = -1;
};

static vector<Recycled_Blk> GlobalHisNodeContainer;
 
class spinlock{
    public:
        void lock() {
            while(lck.test_and_set(std::memory_order_acquire)){}
        };
        void unlock() {
            lck.clear(std::memory_order_release);
        };
    private:
        std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

class Mem_Allocator {
    private:
        Bucket* Bucket_blk = NULL;
        unsigned counter;
        HistoryNode* history_block = NULL;
        Recycled_Blk recycled_nodes;
        unsigned node_counter = 0;
        unsigned recycle_nodecnter = 0;
        unsigned max_level = 0;
    public:
        Mem_Allocator();
        Bucket* Get_bucket();
        HistoryNode* retrieve_his_node();
        unsigned long Request_NodeMem();
        bool Check_Reserve();
};

class Hash_Map {
    private:
        atomic<unsigned long> cur_size;
        int insert_count;
        unsigned long total_size = 0;
        unsigned long max_size = 0;
        unsigned size = 0;
        unsigned thread_cnt = 0;
        size_t node_size = 0;
        unsigned maximum_usage_cnt = 0;
        vector<Bucket*> History_table;
        vector<spinlock> hash_lock;
        vector<Mem_Allocator> Mem_Manager;
    public:
        atomic<long> num_of_waits;
        Hash_Map(size_t size);
        bool isReserve_Empty();
        size_t get_max_size();
        size_t get_cur_size();
        unsigned long get_free_mem();
        uint32_t hash_func(pair<vector<bool>,int> item);
        void calculate_SD();
        void adjust_max_size(int node_count, int GPQ_size);
        void increase_size(size_t size_incre);
        void set_node_t(int node_count, int sub_size, int thread_total);
        void set_up_mem(int thread_num,int node_count);
        void set_sample_period(unsigned runtime_limit);
        void average_size();
        void test_check();
        bool key_compare(const vector<bool>& src, const vector<bool>& dest);
        void analyze_table();
        void print_curmem();
        HistoryNode* insert(pair<boost::dynamic_bitset<>,int>& item,int prefix_cost,int lower_bound, unsigned thread_id, bool backtracked, unsigned depth);
        void table_lock(size_t key);
        void table_unlock(size_t key);
        HistoryNode* retrieve(pair<boost::dynamic_bitset<>,int>& item,int key);
};

#endif