#include "history.hpp"
#include <limits>
#include <vector>

HistoryNode::HistoryNode(int prefix_cost, int bound){
	Entry = {prefix_cost,bound};
}

HistoryNode::HistoryNode(){
    Entry = {-1,-1};
}
