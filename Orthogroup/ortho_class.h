// Created by zhangshengbin on 2021/7/20.
#ifndef _ORTHO_H
#define _ORTHO_H

#include <fstream>
#include <cstring>
#include <iostream>
#include <sys/mman.h>
#include <mutex>
#include <thread>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include "other_function.h"
#include "sliding_mode.h"

using namespace std;
typedef map<string, vector<string> >::iterator my_iter;
static mutex mtx;

class OthoHGT{
public:
    // 构造函数
    OthoHGT(string &file_name, Parse_Tree &obj_tree, bool sd_mode, double min_th, sliding_mode &sl_matrix, map<int, unordered_set<string>> &time_node_map):_ortho_file_name(file_name), _sd_mode(sd_mode), _min_th(min_th), _obj_tree(obj_tree), _sl_matrix(sl_matrix), _time_node_map(time_node_map){}

    // 读取并解析两种类型orthogroups文件
    inline void read_file();

    // 水平基因转移多线程启动函数，建立线程池，启动任务
    inline void statistic_HGT_genes_paraller(int nthreads);


    map<string, vector<string> > ortho_prot;
    map<string, unordered_set<string>> ortho_node_present;
    map<int, set<string>> time_ortho_map;

    // 析构函数
    ~OthoHGT(){}

private:
    string &_ortho_file_name;

    map<int, unordered_set<string>> _time_node_map;

    // 解析的树文件产生的内容
    Parse_Tree &_obj_tree;

    // 滑动模型矩阵
    sliding_mode &_sl_matrix;

    // 文件行数的位数，用来根据行数给各个 Ortho 命名
    int _file_line_figure;

    // 用这个值来记录 orthogroups 文件是哪个软件产生的
    bool _orthofinder2;

    // 记录滑动模式是否开启
    bool _sd_mode;
    // 记录最小阈值
    double _min_th;


    //这两个变量用于统计 HGT 推断的进度
    double _progress = 0;
    double _f_line = 0;


    unordered_map<string ,unordered_map<string, int>> _ortho_genome_num;
    map<string, vector<string> > _ortho_genome_map;

    // 返回文件的行数
    inline void _get_file_line();
    inline void _parse_file(ifstream &infile, vector<string> &, vector<string> &, int, string, string);
    // 水平基因转移推断函数，下面的几个函数都是用于配合它的
    inline void _statistic_HGT_genes(my_iter begin, my_iter end);
    // 每个 ortho 在树上分布深度查找函数
    inline void _ortho_depth_search(string &every_ortho);
    // 这个函数对于任何一个同源家族，从根节点向叶子节点查找其分布大于75%的树
    // 已经查到大于75%的树不再查其子树，转而查其兄弟树
    // 小于75%的树会一路向下查，直到查到叶节点(单个叶节点如果含有为100%，也考虑在内了)为止
    // 叶子节点之间的水平基因转移也被包含在内
    void _leaf_mapping_search(string &tree_name, vector<string> &ortho_genome_vec, vector<string> &satisfied_tree, vector<string> &inspect_vec);
    void _find_all_child(string &ortho_name, string &tree_name);

};

inline void OthoHGT::_get_file_line(){
    size_t length;
    auto f = map_file(_ortho_file_name.c_str(), length);
    auto l = f + length;
    auto f_copy = f;
    while(f && f!=l)
        if((f = static_cast<char*>(memchr(f, '\n', l-f))))
            _f_line++, f++;
    munmap(f_copy, length);
    cout << "We detect there are " << _f_line - 1 << " orthogroups in your orthogroups file." << endl;
    _file_line_figure = give_num_figure(_f_line);
}

inline void OthoHGT::read_file(){
    // 各种变量定义
    ifstream infile(_ortho_file_name.c_str());
    string f_line;

    // 文件没打开，报错并且退出程序
    if(! infile)
        handle_error("Error: Input orthogroups file doesn't find!");

        // 如果文件成功打开，进行下面这些处理
    else{
        // 这一部分表示文件是成功打开的
        cout << "Successful open orthogroups file!" << endl;
        _get_file_line();
        getline(infile, f_line);
        vector <string> f_line_vec = split(f_line, "\t");
        vector <string>::iterator begin_pos = f_line_vec.begin();

        // 下面 orthogroups 文件第一行的第一个 cell 来区分是哪个软件
        // 之后根据软件的不同对文件进行解析

        // 如果由 Orthofinder2 提供
        if(f_line_vec[0] == "Orthogroup"){
            _orthofinder2 = true;
            cout << "We detect the orthogroups file you provide is produced by Orthofinder2." << endl;
            _parse_file(infile, f_line_vec, _obj_tree.tree_leaf_map["WHOLE"], 1, ", ", "");
        }

            // 如果由 proteinortho 提供和 Orthofinder2 处理方式一样，但是一些判断条件有变化
        else{
            _orthofinder2 = false;
            cout << "We detect the orthogroups file you provide is produced by proteinortho." << endl;
            _parse_file(infile, f_line_vec, _obj_tree.tree_leaf_map["WHOLE"], 3,  ",", "*");
        }
        infile.close();
    }
}

inline void OthoHGT::_parse_file(ifstream &infile, vector<string> &fl_vec, vector<string> &leaf_vec, int b_pos, string sep, string empty_f){

    int row_pos = 1, col_pos = 0, name_m_f = 0;
    map<int, string> pos_name_map;
    string line, ori_name;
    //对第一行 f_line 进行处理
    for(auto ex_cell = fl_vec.begin() + b_pos; ex_cell < fl_vec.end(); ex_cell++){
        *ex_cell = strip(*ex_cell);
        ori_name = *ex_cell;
        // 引用，直接把名字改了
        change_variable_name(*ex_cell);

        if(find(leaf_vec.begin(), leaf_vec.end(), *ex_cell) == leaf_vec.end()){
            if(name_m_f == 0){
                cerr << "These names in the orthogroups file can't be find in the tree file:" << endl;
                cerr << ori_name << endl;
                name_m_f = 1;
            }
            else
                cerr << ori_name << endl;
        }
        // 储存改过之后的名字，这个地方考虑了不同软件的起始位点偏移
        // 下面从map里提取对应的基因名称的时候也是用 0 开始
        pos_name_map[col_pos] = *ex_cell;
        col_pos++;
    }


    if(name_m_f == 1)
        exit(255);

    cout << "We detect there are " << col_pos << " genomes in your orthogroups file, and all of them are consistent with the names in species tree." << endl;
    // 下面是对第一行之后的内容进行处理
    while(getline(infile, line)){
        // 第一行跳过去前面处理过了
        col_pos = 0;
        string OG_name = give_OG_name(row_pos, _file_line_figure);
        vector<string> line_vec = split(line, "\t");
        for(auto ex_cell = (line_vec.begin() + b_pos); ex_cell < line_vec.end(); ex_cell++) {
            if (*ex_cell != "*" and ex_cell->size() != 0) {
                *ex_cell = strip(*ex_cell);
                _ortho_genome_map[OG_name].push_back(pos_name_map[col_pos]);
                vector<string> prot_vec = split(*ex_cell, sep);
                _ortho_genome_num[OG_name].insert({pos_name_map[col_pos], prot_vec.size()});
                for (auto &i:prot_vec) {
                    strip(i);
                    change_variable_name(i);
                    ortho_prot[OG_name].push_back(i);
                }
            }
            col_pos++;
        }
        row_pos++;
    }
}

inline void OthoHGT::_ortho_depth_search(string &every_ortho){
    // satisfied_tree, inspect_vec 两个 vec 都是每个 ortho 更新一次
    vector<string> ortho_genome_vec = _ortho_genome_map[every_ortho], satisfied_tree, inspect_vec;
    double total_leaf = _obj_tree.tree_leaf_map["WHOLE"].size(), num_genome = ortho_genome_vec.size();

    // 如果这个 ortho 含有基因组数量太少，先找到这些基因组的 LCA，再从 LCA 往下找

    if(num_genome/total_leaf > 0.01){
        string temp = "WHOLE";
        _leaf_mapping_search(temp, ortho_genome_vec, satisfied_tree, inspect_vec);
    }
    else{
        string ancestor = _obj_tree.find_LCA(ortho_genome_vec);
        _leaf_mapping_search(ancestor, ortho_genome_vec, satisfied_tree, inspect_vec);
    }


    for (auto &i:satisfied_tree) {
        if (_obj_tree.check_leaf_node(i)) {
            mtx.lock();
            for(auto &j:_time_node_map){
                if(j.second.find(i) != j.second.end())
                    time_ortho_map[j.first].insert(every_ortho);
            }
            ortho_node_present[every_ortho].insert(i);
            mtx.unlock();
        }
        else {
            mtx.lock();
            for(auto &j:_time_node_map) {
                if (j.second.find(i) != j.second.end())
                    time_ortho_map[j.first].insert(every_ortho);
            }
            ortho_node_present[every_ortho].insert(i);
            mtx.unlock();
            _find_all_child(every_ortho, i);
        }
    }
}

inline void OthoHGT::statistic_HGT_genes_paraller(int nthreads){
    int max_threads = thread::hardware_concurrency();
    // 构建一个线程池
    vector <thread> threads(nthreads);
    map<string, vector<string> > order_change_map;
    // 计算每个线程池里有多少个任务
    int grainsize = _ortho_genome_map.size()/nthreads + 1;
    cout << setiosflags(ios::fixed) << setprecision(2);

    if(nthreads > max_threads)
        cerr << "Warining: The number of threads you requested exceeds the maximum currently available on your computer, "
                "which may cause the program to slow down" << endl;

    move_map_order(_ortho_genome_map, order_change_map, nthreads, grainsize);
    // 对线程池进行任务分配
    map<string, vector<string> >::iterator temp_iter = order_change_map.begin();
    auto work_iter_b = temp_iter;
    for(int i = 0; i < grainsize; i++)
        temp_iter++;
    auto work_iter_e = temp_iter;


    for(auto it = threads.begin(); it != threads.end() - 1; ++it){
        // 原理就是把 map 中的数据拆开了，一部分一部分进行处理
        // 等号右边创建了一个thread对象，将其拷贝初始化给线程池中的某个thread
        *it = thread(&OthoHGT::_statistic_HGT_genes, this, work_iter_b, work_iter_e);
        auto temp_iter2 = work_iter_b;
        work_iter_b = work_iter_e;
        for(int i = 0; i < grainsize; i++)
            temp_iter2++;
        work_iter_e = temp_iter2;
    }
    // 最后一个线程多处理一部分任务，把余数也处理了
    threads.back() = thread(&OthoHGT::_statistic_HGT_genes, this, work_iter_b, order_change_map.end());

    // 等待各个线程运行完毕
    for(auto &i:threads)
        i.join();

}

inline void OthoHGT::_statistic_HGT_genes(my_iter begin, my_iter end){
    for(auto it = begin; it != end; ++it) {
        string ortho_name = it->first;
        mtx.lock();
        _progress++;
        cout << '\r' << "now: " << _progress / _f_line;
        mtx.unlock();
        _ortho_depth_search(ortho_name);
    }
}

#endif