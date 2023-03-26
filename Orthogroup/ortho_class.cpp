// Created by zhangshengbin on 2021/7/20.

#include <algorithm>
#include "ortho_class.h"

void  OthoHGT::_find_all_child(string &ortho_name, string &tree_name){
    for(auto &i:_obj_tree.parent_child_map[tree_name]){
        mtx.lock();
        for(auto &j:_time_node_map) {
            if (j.second.find(i) != j.second.end())
                time_ortho_map[j.first].insert(ortho_name);
        }
        ortho_node_present[ortho_name].insert(i);
        mtx.unlock();
        if(!_obj_tree.check_leaf_node(i))
            _find_all_child(ortho_name, i);
    }
}


void OthoHGT::_leaf_mapping_search(string &tree_name, vector<string> &ortho_genome_vec, vector<string> &satisfied_tree, vector<string> &inspect_vec){
    double tree_map_n = 0.0, total_leaf_n;
    inspect_vec.push_back(tree_name);
    double threshold = _min_th;
    if(_obj_tree.check_leaf_node(tree_name)){
        total_leaf_n = 1.0;
        if(find(ortho_genome_vec.begin(), ortho_genome_vec.end(), _obj_tree.leaf_code_to_name_map[tree_name]) != ortho_genome_vec.end())
            tree_map_n = 1.0;
    }
    else{
        total_leaf_n = _obj_tree.tree_leaf_map[tree_name].size();
        for(auto &eve_leaf:_obj_tree.tree_leaf_map[tree_name])
            if(find(ortho_genome_vec.begin(), ortho_genome_vec.end(), eve_leaf) != ortho_genome_vec.end())
                tree_map_n++;
    }
    // 如果打开了滑动模式，阈值每次都更新，否则就是不变的了
    if(_sd_mode){
        threshold = _sl_matrix.node_threshold_map[tree_name];
    }
    // 如果一棵树完全不含有目标同源家族，或者是含有基因组数量满足要求，那么只查它的兄弟树
    if(tree_map_n/total_leaf_n >= threshold || tree_map_n == 0){
        if(tree_map_n/total_leaf_n >= threshold)
            satisfied_tree.push_back(tree_name);
        if(tree_name != "WHOLE"){
            for(auto &brother_tree:_obj_tree.parent_child_map[_obj_tree.child_parent_map[tree_name]])
                if(find(inspect_vec.begin(), inspect_vec.end(), brother_tree) == inspect_vec.end())
                    _leaf_mapping_search(brother_tree, ortho_genome_vec, satisfied_tree, inspect_vec);
        }
    }
    // 如果不是上面的情况，一直往下查，查它的子树
    else{
        if(!_obj_tree.check_leaf_node(tree_name))
            for(auto &child_tree:_obj_tree.parent_child_map[tree_name])
                if(find(inspect_vec.begin(), inspect_vec.end(), child_tree) == inspect_vec.end())
                    _leaf_mapping_search(child_tree, ortho_genome_vec, satisfied_tree, inspect_vec);
    }
}
