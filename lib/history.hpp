#ifndef HISTORY_H
#define HISTORY_H

#include <iostream>
#include <vector>

using namespace std;

struct HistoryContent {
	int prefix_cost = -1;
	int lower_bound = -1;
};

class HistoryNode {
	public:
		unsigned usage_cnt = 0;
		HistoryContent Entry;
        HistoryNode(int prefix_cost, int bound);
		HistoryNode();
};

#endif