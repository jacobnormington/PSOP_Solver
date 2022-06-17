#ifndef HISTORY_H
#define HISTORY_H

#include <iostream>
#include <vector>
#include <atomic>
#include <cstdint>
#include <mutex>

using namespace std;

class HistoryNode;

struct HistoryContent {
	int32_t prefix_cost = -1;
	int32_t lower_bound = -1;
};

struct Active_Node {
    Active_Node* parent = NULL;
    HistoryNode* history_link = NULL;
	int32_t total_children_cnt = 0;
    int32_t cur_children_cnt = 0;
    atomic<int32_t> cur_threadID;
	atomic<bool> deprecated;
	mutex node_lock;
};

struct HistoryNode {
	atomic<bool> explored;
	atomic<uint8_t> active_threadID;
	atomic<HistoryContent> Entry;
};

#endif