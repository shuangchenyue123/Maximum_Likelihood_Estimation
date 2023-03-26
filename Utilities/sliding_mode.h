// Created by zhangshengbin on 2021/11/26
// 为刘婕献上鲜血
#ifndef _SL_H
#define _SL_H

#include "tree_class.h"
#include <cmath>
#include <unordered_map>

using namespace std;

// 找到每个叶子节点在树上的深度（这里指从根节点到此叶节点有多少代）

class sliding_mode{
public:
    sliding_mode(Parse_Tree &tree, double min_th):_obj_tree(tree), _min_th(min_th){};

    // 计算树上各个节点的阈值，保存在 node_threshold_map 字典中
    inline void calculate_every_node_threshold();

    unordered_map<string, double> node_threshold_map;

    // 保存计算得到的基因丢失速率
    double loss_p = -1;

private:

    // 这个函数来找整棵树所有父子节点距离的平均值
    // 当然 不 包括WHOLE到root，原因在下面函数定义之前写了
    inline void _find_average_dis(string &node_name);

    // 下面这两个函数用来生成所有的系数矩阵
    inline void _make_leaf_loss_possibility_matrix_auxiliary(string &node_name, string &leaf_name, string &target_node);
    inline void _make_leaf_loss_possibility_matrix(string &node_name);

    // 解方程，求出平均基因丢失速率
    inline void _solve_equations();


    map<string, vector<double>> _leaf_k_map;
    Parse_Tree &_obj_tree;
    double _total_dis = 0, _child_num = 0, _ave_dis;
    double _min_th;
    string temp_root = "WHOLE";

};

// >=3 带根形式, 忽略了WHOLE上面的那个root
// 不需要看那个，因为从比例来说root和WHOLE没有区别，所以HGT推断皆是从WHOLE开始而非root，所以就更不需要了
// 我们求平均值，不加这一层也一样求，这样各个形式的树可以统一操作，而且 >=3带根 这种树其实很少很少
inline void sliding_mode::_find_average_dis(string &node_name) {
    for(auto &i:_obj_tree.parent_child_map[node_name]){
        // 不是叶子节点，距离加和，节点数加1，递归
        if(!_obj_tree.check_leaf_node(i)) {
            _total_dis += _obj_tree.give_tree_evolution_distance_father(node_name);
            _child_num ++;
            _find_average_dis(i);
        }
        // 是叶子节点，距离加和，节点数加1
        else {
            _total_dis += _obj_tree.give_tree_evolution_distance_father(node_name);
            _child_num ++;
        }
    }
}

inline void sliding_mode::_make_leaf_loss_possibility_matrix_auxiliary(string &node_name, string &leaf_name, string &target_node){
    // 找一下父节点是谁，作为递归参数
    string father = _obj_tree.child_parent_map[node_name];

    // 计算父子距离并且和平均距离做比较，生成系数，这里的系数无所谓有是这个节点还是其父节点的系数的问题
    // 因为对同一个叶子节点递归，所有系数生成一个向量
    // 不同的叶子节点的向量组成一个矩阵
    // 对矩阵进行运算得到的是递归终点的阈值，所以这里不需要区分
    double dis = _obj_tree.give_tree_evolution_distance_father(node_name);
    double per = dis/_ave_dis;

    // 控制系数范围，我们认为基因丢失速度的最大差值是100倍，如果上面模型超出这个范围，人工调整回来
    if(per > 10)
        per = 9.9999;
    else if(per < 0.1)
        per = 0.1001;
    // leaf_name这里全都是用的代号，不是名字
    _leaf_k_map[leaf_name].push_back(per);

    if(father != target_node)
        _make_leaf_loss_possibility_matrix_auxiliary(father, leaf_name, target_node);
}

// 系数矩阵的生成算法
inline void sliding_mode::_make_leaf_loss_possibility_matrix(string &node_name){
    // 根据另一个函数的结果，先计算一个平均距离
    _ave_dis = _total_dis/_child_num;

    // 找到我们要计算节点对应的所有叶子节点，遍历
    for(auto leaf:_obj_tree.tree_leaf_map[node_name]){
        // 这里要把叶子节点的名字切换成代号，辅助函数里的字典只收代号
        leaf = _obj_tree.leaf_name_to_code_map[leaf];
        // 传递的三个参数用途依次是递归开始的节点，用于作储存字典的键值，用于作递归终点
        _make_leaf_loss_possibility_matrix_auxiliary(leaf, leaf, node_name);
    }
}


inline void sliding_mode::_solve_equations() {
    // 解方程，通过简化运算可知当i在[0.001,0.1]之间的时候，整个矩阵计算的结果大约是[0.28,0.99]
    // 我们以0.0005作为增量，当结果位于[min_th-0.01, min_th+0.01]时，结束运算
    for(double i = 0.001; i <= 0.1; i += 0.0005) {
        // 各个向量的和，也是最终结果
        double temp_total = 0;
        for (auto leaf:_leaf_k_map) {
            // 每个向量都是一个乘积运算
            double temp_leaf = 1;
            //顺带把所有leaf的阈值都置换成1，这个和解方程没关系，这里捎带着了
            node_threshold_map[leaf.first] = 1;
            for (auto &per:leaf.second) {
                double temp = 1 - per * i;
                temp_leaf *= temp;
            }
            temp_total += temp_leaf;
        }
        temp_total = temp_total/_leaf_k_map.size();
        if(temp_total < _min_th + 0.01 and temp_total > _min_th - 0.01){
            loss_p = i;
            _leaf_k_map.clear();
            break;
        }
    }
}


inline void sliding_mode::calculate_every_node_threshold(){
    // 先把标准丢失速率解出来
    _find_average_dis(temp_root);
    _make_leaf_loss_possibility_matrix(temp_root);
    _solve_equations();
    if(loss_p != -1)
        cout << "The standard gene loss possibility we detect is: " << loss_p << "." << endl;
    else
        handle_error("The min_threshold you provide is too low");

    // 对树上的每一个节点求阈值，根节点不用求，是用户提供的
    // 叶子节点一定都是1，我们上面已经设置过了
    for(auto &node:_obj_tree.parent_child_map){
        // 如果是根节点我们直接跳过
        if(node.first == "WHOLE")
            continue;

        // 除了根节点以外，每一个都计算一个阈值
        // 针对每一个节点，先生成一个系数矩阵
        _make_leaf_loss_possibility_matrix(const_cast<string &>(node.first));

        double temp_total = 0;
        // 根据标准丢失速率和系数矩阵进行运算
        for (auto &leaf:_leaf_k_map) {
            double temp_leaf = 1;
            for (auto &per:leaf.second) {
                double temp = 1 - per * loss_p;
                temp_leaf *= temp;
            }
            temp_total += temp_leaf;
        }
        temp_total = temp_total/_leaf_k_map.size();
        node_threshold_map[node.first] = temp_total;
        _leaf_k_map.clear();
    }

    // 上面那个没有计算叶子节点，这里要补上，叶子节点实际上全部都是1
    for(auto &leaf:_obj_tree.tree_leaf_map["WHOLE"]) {
        node_threshold_map[_obj_tree.leaf_name_to_code_map[leaf]] = 1;
    }

    node_threshold_map["WHOLE"] = _min_th;
}


#endif