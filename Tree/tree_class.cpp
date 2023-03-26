// Created by zhangshengbin on 2021/7/20.


#include <bits/stdc++.h>
#include <algorithm>
#include "tree_class.h"

void Parse_Tree::parse_tree(string &tree_text, string &tree_name){
    string::iterator every_char;
    int left_brackets_num = -1;
    int right_brackets_num = 0;
    string child_tree_text;
    int surplus_bracket_num = _statistic_bracket_num(tree_text);

    tree_text = tree_text.substr(surplus_bracket_num, tree_text.size()-surplus_bracket_num);
    if(tree_name == _root)
        name_text_map[tree_name] = tree_text;
    for(every_char = tree_text.begin(); every_char < tree_text.end(); every_char++){
        switch (*every_char) {
            case '(':
                left_brackets_num += 1;
                break;
            case ')':
                right_brackets_num += 1;
                if(right_brackets_num > left_brackets_num) {
                    string child_tree_name = _make_relationship(child_tree_text, tree_name);
                    parse_tree(child_tree_text, child_tree_name);
                }
                break;
            case ',':
                if(left_brackets_num == right_brackets_num) {
                    string child_tree_name = _make_relationship(child_tree_text, tree_name);
                    parse_tree(child_tree_text, child_tree_name);
                    child_tree_text = "";
                    continue;
                }
                else
                    break;
            case ' ': case '\t': case '\n':
                continue;
            default:
                break;
        }
        child_tree_text += *every_char;
    }
}

void Parse_Tree::find_tree_leaf(string &tree_name){
    for(auto &it:parent_child_map[tree_name]){
        if(check_leaf_node(it)){
            string tree_text = name_text_map[it];
            string temp;
            if(_time_iq)
                temp = tree_text.substr(0,tree_text.find_last_of("["));
            else
                temp = tree_text.substr(0,tree_text.find_last_of(":"));
            change_variable_name(temp);
            leaf_name_to_code_map[temp] = it;
            leaf_code_to_name_map[it] = temp;
            tree_leaf_map[tree_name].push_back(temp);
        }
        else{
            find_tree_leaf(it);
            if(tree_leaf_map.count(tree_name) != 0)
                tree_leaf_map[tree_name].insert(tree_leaf_map[tree_name].end(), tree_leaf_map[it].begin(), tree_leaf_map[it].end());
            else
                tree_leaf_map[tree_name] = tree_leaf_map[it];
        }
    }
}

void Parse_Tree::find_every_tree_child(string &tree_name, map<string, bool> &temp_map){
    for(auto &child_tree:parent_child_map[tree_name]){
        temp_map[child_tree] = true;
        if(! check_leaf_node(child_tree))
            find_every_tree_child(child_tree, temp_map);
    }
}

// 这个函数的内存管理有问题，一会儿回来改
string Parse_Tree::find_LCA(vector<string> &leaf_group, int search_t){
    vector<string> father_list, father_list_clean;
    map<string, bool> temp_map;
    if(leaf_group.size() == 1)
        return leaf_group[0];
    else{
        // 遍历每个叶子节点
        for(auto &it_lf:leaf_group){
            string leaf_code, parent_code;
            // 第一次查找的时候要先找到叶子节点对应的代号才能进入 child_parent_map 字典
            if(search_t == 1){
                leaf_code = leaf_name_to_code_map[it_lf];
                parent_code = child_parent_map[leaf_code];
            }
            else
                parent_code = child_parent_map[it_lf];

            // 把本次每一个待查找节点的父节点放进一个向量里
            // 放进去之前检查一下父节点之间有没有重复的，重复的就不放进去了
            if(find(father_list.begin(), father_list.end(), parent_code) == father_list.end())
                father_list.push_back(parent_code);
            // 查找本次所有父节点对应的后代节点
            find_every_tree_child(parent_code, temp_map);
        }
        // 继续检查之前的 father_list 里，有没有 father 其实是其它 father 的child
        // 这种的被剔除掉，不进行下一轮的查验
        for(auto &it:father_list)
            if(temp_map.count(it) == 0)
                father_list_clean.push_back(it);
        // 此处只区分是不是 1，后面都用 2 表示就可以了
        return find_LCA(father_list_clean, 2);
    }
}

// 计算每一个点到根节点的距离
// 除了 WHOLE 没有，该有的都有，WHOLE 不需要有
// >=3 带根形式，区别是这个形式的 WHOLE 能取到 sh/bootstrap/evol_dis
// 因为这种形式大家到根节点距离一样，因此可以舍掉，二叉树的根节点是可以被计算在内的
// 这样 >=3 无根的形式不存在计算上的问题，但是由于根不确定，导致大家到 WHOLE 进化距离有变数！
// 我们这里不做多余处理，默认根在哪个点就是哪个点了，反正也给出了警告信息

// 计算各棵树的检验值
// WHOLE大家都一样，就不用算了
// 叶子节点的返回值问题已经在对应的 self 函数中考虑过了
void Parse_Tree::make_tree_attribution_map(string &tree_name, double total_dis, double total_test, int layers){
    for(auto &tree:parent_child_map[tree_name]){
        double total_dis_copy = total_dis, total_test_copy = total_test;
        double temp_time = give_tree_evolution_distance_father(tree);
        double temp_test = give_test_result_self(tree);
        total_dis_copy += temp_time;
        // 直接把浮点数的误差去掉
        total_test_copy += temp_test;
        layers ++;
        tree_evolution_map[tree] = total_dis_copy;
        tree_test_map[tree] = total_test_copy/layers;
        if(!check_leaf_node(tree))
            make_tree_attribution_map(tree, total_dis_copy, total_test_copy, layers);
    }
}

