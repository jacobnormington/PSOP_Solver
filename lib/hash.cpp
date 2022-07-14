#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <sys/sysinfo.h>
#include <boost/dynamic_bitset.hpp>
#include "hash.hpp"

Mem_Allocator::Mem_Allocator() {
    Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
    history_block = new HistoryNode[HIS_BLK_SIZE];
    counter = 0;
    node_counter = 0;
}

HistoryNode* Mem_Allocator::retrieve_his_node() {
    if (node_counter >= HIS_BLK_SIZE || history_block == NULL) {
        history_block = new HistoryNode[HIS_BLK_SIZE];
        node_counter = 0;
    }
    HistoryNode* node = history_block + node_counter;
    node_counter++;
    return node;
}


Bucket* Mem_Allocator::Get_bucket() {
    if (counter == BUCKET_BLK_SIZE || Bucket_blk == NULL) {
        Bucket_blk = new Bucket[BUCKET_BLK_SIZE];
        counter = 0;
    }
    Bucket* bucket = Bucket_blk + counter;
    counter++;
    return bucket;
}

Hash_Map::Hash_Map(int size) {
    this->size = size;
    struct sysinfo info;
	if(sysinfo(&info) != 0){
        cout << "can't retrieve sys mem info\n";
		exit(1);
	}
    max_size = ((double)info.freeram * MEMORY_RESTRIC) - (size * 4);
    cur_size = 0;
    History_table.resize(size);
    
    for (int i = 0; i < size; i++) History_table[i] = NULL;
}

size_t Hash_Map::get_max_size() {
    return max_size;
}

size_t Hash_Map::get_cur_size() {
    return cur_size;
}

void Hash_Map::increase_size(size_t size_incre) {
    cur_size += size_incre;
    return;
}

void Hash_Map::set_node_t(int node_count) {
    node_size = 2*sizeof(HistoryNode*) + sizeof(HistoryNode) + sizeof(boost::dynamic_bitset<>) + sizeof(int) + (size_t)max(node_count/8,64);
    return;
}

void Hash_Map::average_size() {
    long unsigned max = 0;
    long unsigned sum = 0;
    long unsigned item = 0;
    for (unsigned i = 0; i < size; i++) {
        if (History_table[i] != NULL) {
            if (History_table[i]->size() > 0) item++;
            if (History_table[i]->size() > max) max = History_table[i]->size();
            sum +=  History_table[i]->size();
        }
    }
    cout << "Total non empty entry is " << item << endl;
    cout << "Total history table size is " << sum << endl;
    cout << "Maximum bracket size is " << max << endl;
    cout << "Average bracket size is " << float(sum)/float(item) << endl;
}
    
bool Hash_Map::insert(pair<boost::dynamic_bitset<>,int>& item,int prefix_cost,int lower_bound, int best_cost) {
    size_t val = hash<boost::dynamic_bitset<>>{}(item.first);
    int key = (val + item.second) % size;

    HistoryNode* node = Mem_Manager.retrieve_his_node();
    node->Entry = {prefix_cost,lower_bound};

    if (History_table[key] == NULL) {
        History_table[key] = Mem_Manager.Get_bucket();
        History_table[key]->push_back(make_pair(item,node));
        cur_size += (node_size + sizeof(Bucket));
        return true;
    }

    History_table[key]->push_back(make_pair(item,node));
    cur_size += node_size;
    return true;
}


HistoryNode* Hash_Map::retrieve(pair<boost::dynamic_bitset<>,int>& item,int key) {
    if (History_table[key] == NULL) return NULL;
    else if (History_table[key]->size() == 1) {
        if (item.second == History_table[key]->begin()->first.second && item.first == History_table[key]->begin()->first.first) {
            return History_table[key]->begin()->second;
        }
        else return NULL;
    }

    for (auto iter = History_table[key]->begin(); iter != History_table[key]->end(); iter++) {
        if (item.second == iter->first.second && item.first == iter->first.first) {
            return iter->second;
        }
    }

    return NULL;
}


void Hash_Map::analyze_usage() {
    unsigned long total_entry = 0;
    vector<int> record = vector<int>(maximum_usage_cnt+1,0);
    for (unsigned i = 0; i < size; i++) {
        
        if (History_table[i] != NULL) {
            total_entry += History_table[i]->size();
            for (auto iter = History_table[i]->begin(); iter != History_table[i]->end(); iter++) {
                record[iter->second->usage_cnt]++;
            }
        }   
    }

    cout << "Total number of entries is " << total_entry << endl;

    for (unsigned i = 0; i < maximum_usage_cnt; i++) {
        if (record[i] != 0) cout << i << "," << record[i] << endl;
    }

}