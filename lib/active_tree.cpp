#include "active_tree.hpp"
#include "solver.hpp"
#include <climits>
#include <vector>
#include <queue>
#include <deque>
#include <mutex>

mutex active_map_lock;

Active_Path::Active_Path() {
    size = 0;
    ptr_cnt = 0;
}

Active_Path::Active_Path(int num, int arr_size) {
    size = num;
    path = vector<Active_Node*>(size);
    for (unsigned i = 0; i < path.size(); i++) {
        path[i] = new Active_Node();
        path[i]->cur_children_cnt = -1;
        path[i]->total_children_cnt = -1;
        path[i]->history_link = NULL;
        path[i]->parent = NULL;
    }
    ptr_cnt = arr_size - 2;
}

Active_Node* Active_Path::get_element(int index) {
    if (index < size && index >= 0) return path[index];
    return NULL;
}

Active_Node* Active_Path::back() {
    if (ptr_cnt < size && ptr_cnt >= 0) return path[ptr_cnt];
    return NULL;
}

HistoryNode* Active_Path::back_history() {
    if (ptr_cnt < size && ptr_cnt >= 0) return path[ptr_cnt]->history_link;
    return NULL;
}


size_t Active_Path::get_size() {
    return ptr_cnt;
}

void Active_Path::set_threadID(int id, int total_id) {
    thread_id = id;
    thread_total = total_id;
}

bool Active_Path::push_back(int children_num, HistoryNode* his_node, Active_Allocator* Allocator) {
    ptr_cnt++;
    if (ptr_cnt >= size) {
        ptr_cnt = size - 1;
        cout << "push_back error" << endl;
        exit(1);
    }
    path[ptr_cnt]->total_children_cnt = children_num;
    path[ptr_cnt]->cur_threadID = thread_id;
    path[ptr_cnt]->cur_children_cnt = 0;

    if (his_node != NULL && his_node->active_threadID == thread_id) {
        path[ptr_cnt]->history_link = his_node;
        his_node->explored = false;
    }
    else path[ptr_cnt]->history_link = NULL;

    path[ptr_cnt]->parent = path[ptr_cnt - 1];
    return true;
}

bool Active_Path::pop_back(Active_Allocator* Allocator, bool out_dated) {
    if (ptr_cnt <= 0) {
        ptr_cnt = 0;
        cout << "pop_back error" << endl;
        exit(1);
    }

    path[ptr_cnt]->node_lock.lock();
    if (path[ptr_cnt]->total_children_cnt <= path[ptr_cnt]->cur_children_cnt) {
        if (path[ptr_cnt]->parent != NULL) {
            path[ptr_cnt]->parent->node_lock.lock();
            path[ptr_cnt]->parent->cur_children_cnt++;
            if ((path[ptr_cnt]->parent->total_children_cnt <= path[ptr_cnt]->parent->cur_children_cnt) && path[ptr_cnt]->parent->cur_threadID < 0) {
                path[ptr_cnt]->parent->node_lock.unlock();
                upward_propagation(Allocator);
            }
            else path[ptr_cnt]->parent->node_lock.unlock();
        }
        if (path[ptr_cnt]->history_link != NULL && path[ptr_cnt]->history_link->active_threadID == thread_id) {
            path[ptr_cnt]->history_link->explored = true;
        }
        path[ptr_cnt]->deprecated = false;
        path[ptr_cnt]->node_lock.unlock();
    }
    else {
        path[ptr_cnt]->cur_threadID = -1;
        path[ptr_cnt]->deprecated = out_dated;
        path[ptr_cnt]->node_lock.unlock();
        path[ptr_cnt] = Allocator->assign_node();
    }
    ptr_cnt--;

    return true;
}

bool Active_Path::check_deprecation_status(int level) {
    if (level > size || level <= 0 || path[level] == NULL) return false;
    if (path[level]->deprecated) return true;
    else return false;

    return false;
}

bool Active_Path::incre_children_cnt(Active_Allocator* Allocator) {
    if (ptr_cnt > 0 && ptr_cnt < size) {
        path[ptr_cnt]->node_lock.lock();
        path[ptr_cnt]->cur_children_cnt++;
        if (path[ptr_cnt]->total_children_cnt <= path[ptr_cnt]->cur_children_cnt && path[ptr_cnt]->cur_threadID < 0) {
            if (path[ptr_cnt]->parent != NULL) {
                if (path[ptr_cnt]->parent != path[ptr_cnt-1]) {
                    cout << endl;
                    cout << "ACTIVE_TREE ERROR" << endl;
                    cout << "ptr cnt is " << ptr_cnt << endl;
                    cout << "/////Path[ptr_cnt - 1]/////" << endl;
                    cout << path[ptr_cnt-1]->total_children_cnt << path[ptr_cnt-1]->cur_children_cnt << endl;
                    cout << "/////Path[ptr_cnt].parent/////" << endl;
                    cout << path[ptr_cnt]->parent->total_children_cnt << path[ptr_cnt-1]->parent->cur_children_cnt << endl;
                    cout << endl;

                    enable_threadstop = false;
                }
                path[ptr_cnt]->parent->node_lock.lock();
                path[ptr_cnt]->parent->cur_children_cnt++;
                if (path[ptr_cnt]->parent->total_children_cnt <= path[ptr_cnt]->parent->cur_children_cnt && path[ptr_cnt]->parent->cur_threadID < 0) {
                    path[ptr_cnt]->parent->node_lock.unlock();
                    upward_propagation(Allocator);
                }
                else path[ptr_cnt]->parent->node_lock.unlock();
            }
            
            if (path[ptr_cnt]->history_link != NULL && path[ptr_cnt]->history_link->active_threadID == thread_id) {
                path[ptr_cnt]->history_link->explored = true;
            }
            Allocator->delete_node(path[ptr_cnt]);
        }
        path[ptr_cnt]->node_lock.unlock();
    }
    else return false;

    return false;
}

bool Active_Path::assign_hislink(HistoryNode* link) {
    if (path[ptr_cnt] != NULL) path[ptr_cnt]->history_link = link;
    return false;
}

int Active_Path::upward_propagation(Active_Allocator* Allocator) {
    //cout << "upward_propagation initialized" << endl;
    
    for (int i = ptr_cnt - 1; i > 0; i--) {
        path[i]->node_lock.lock();
        if (path[i]->total_children_cnt != -1 && path[i]->total_children_cnt <= path[i]->cur_children_cnt && path[i]->cur_threadID < 0) {
            if (path[i]->parent != NULL) {
                path[i]->parent->node_lock.lock();
                path[i]->parent->cur_children_cnt++;
                path[i]->parent->node_lock.unlock();
            }
            else {
                //print_curpath();
                path[i]->node_lock.unlock();
                return i;
            }
    
            if (path[i]->history_link != NULL) {
                path[i]->history_link->explored = true;
                path[i]->history_link = NULL;
            }
            Allocator->delete_node(path[i]);
        }
        else {
            path[i]->node_lock.unlock();
            return i;
        }
        path[i]->node_lock.unlock();
    }

    return 1;
}

void Active_Path::reconstruct_path(Active_Path& partial_path, Active_Allocator* Allocator) {
    //cout << "reconstructed...." << endl;
    ptr_cnt = partial_path.get_size();

    for (int i = 0; i <= ptr_cnt; i++) {
        /*
        cout << "[" << partial_path.get_element(i)->total_children_cnt << "," 
                     << partial_path.get_element(i)->cur_children_cnt <<  "," 
                     << partial_path.get_element(i)->cur_threadID << "],";
        */
        path[i] = partial_path.get_element(i);
        if (i > 0) path[i]->parent = path[i-1];
    }

    for (int i = ptr_cnt + 1; i < size; i++) {
        path[i] = Allocator->assign_node();
        path[i]->cur_children_cnt = -1;
        path[i]->total_children_cnt = -1;
        path[i]->cur_threadID = -1;
        path[i]->history_link = NULL;
        path[i]->parent = NULL;
        path[i]->deprecated = false;
    }
    //cout << endl;
    return;
}

void Active_Path::clear(Active_Allocator* Allocator) {
    for (int i = size-1; i > ptr_cnt; i--) {
        Allocator->delete_node(path[i]);
    }
    //cout << "ptr_cnt = " << ptr_cnt << endl;
    return;
}

void Active_Path::print_curpath() {
    for (int i = 0; i < size; i++) cout << "[" << path[i]->total_children_cnt 
                                            << "," << path[i]->cur_children_cnt 
                                            << "," << path[i]->cur_threadID << "]" << ",";
    cout << endl;
    
    return;
}