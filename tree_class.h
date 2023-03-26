// Created by zhangshengbin on 2021/7/20.
// 我以血肉供养挚爱
// 为刘婕献上鲜血
#ifndef _TREE_CLASS_H
#define _TREE_CLASS_H

#include <regex>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include "other_function.h"


using namespace std;

// 下面这个类用于解析树文件
class Parse_Tree{
public:
    Parse_Tree(string &file_name, bool bootstrap, bool sh, bool time_iq):_tree_file_name(file_name), t_bootstrap(bootstrap), t_sh(sh), _time_iq(time_iq){}

    // 声明读取树文件的函数
    inline void read_tree();

    // 整个用来对树进行解析，是这个类的核心
    void parse_tree(string &tree_text, string &tree_name);

    // 下面是一些功能函数，针对解析过的树，实现一些目标功能

    // 下面这个函数专门用于判断树的模式，返回值设置为1、2、3
    // 1 代表 >=3 不带根
    // 2 代表 >=3 带根
    // 3 代表 <3 带根
    inline void judge_tree_root();

    // 用来检查某一棵树是否是叶子节点
    inline bool check_leaf_node(string &tree_name);

    // 这个函数给出每颗树到其父节点的距离，只在下一个函数中套用，不单独使用
    // 对于 treePL 来说，时间=距离
    // 由于下一个函数的特性，tree_name 永远不会是 WHOLE
    inline double give_tree_evolution_distance_father(string &tree_name);
    // 统计sh检验和bootstrap检验值
    inline double give_test_result_self(string &tree_name);
    // 这个函数给出任一一对祖先和后代节点之间的进化距离，树定根的问题在函数定义部分有注释进行了讨论
    void make_tree_attribution_map(string &tree_name, double total_dis = 0, double total_test = 0, int layers = 0);


    // 找到各个树的代号和其对应的叶子节点名称，找到叶子节点代号和其名称的对应关系
    void find_tree_leaf(string &tree_name);

    // 针对iqtree2时间树，单独设置一个函数查找其时间
    inline double search_iqtree2_time(string &tree_name);

    // 用来查找某棵树上某个节点全部的子节点，一直到叶子节点, 主要用来配合下面这个函数
    void find_every_tree_child(string &tree_name, map<string, bool> &temp_map);
    // 用来查找树上某些叶子节点最近的共同祖先
    string find_LCA(vector<string> &leaf_group, int search_t = 1);

    // 这两个函数用于程序测试，没有其它用途
    inline void print_tree_leaf();
    inline void print_child_tree();

    bool t_bootstrap;
    bool t_sh;
    double oat;
    unordered_map<string, string> name_text_map;
    map<string, vector<string> > parent_child_map;
    map<string, string> child_parent_map;
    map<string, vector<string> > tree_leaf_map;
    unordered_map<string, string> leaf_name_to_code_map;
    unordered_map<string, string> leaf_code_to_name_map;
    unordered_map<string, double> tree_evolution_map;
    unordered_map<string, double> tree_test_map;

    // 没有 new/malloc 的时候析构函数不需要手动进行处理
    ~Parse_Tree(){}

    // 用的是内部变量的形式，但是为了外部第一次调用这里还是给到一个外部变量
    string _tree_text;

private:
    const string _tree_file_name;
    const string _root = "WHOLE";
    bool _time_iq;
    const double _scale = 0.000001;
    regex _time_pattern{R"(&date=\"(-?[\.\d]*)\"\]:[\d\.]*$)"};
    regex _test_double_pattern{R"(([\d\.]*)\/([\d\.]*)[^\)]*$)"};
    regex _test_single_pattern{R"(([\d\.]*)[^\)]*$)"};

    // 下面的三个函数是用来配合上面的树解析函数来设立的
    // 随机生成三个不重复的大写字母作为各个树的名字
    inline string _make_tree_name(string &tree_name);
    inline string _make_relationship(string &tree_text, string &father_tree_name);
    inline int _statistic_bracket_num(string &tree_text);
};

// 可以同时处理 newick 和 nexus 两种格式的文件
inline void Parse_Tree::read_tree(){
    string line;
    // 打开输入文件
    bool store_flag = 0;
    ifstream infile(_tree_file_name);
    if(! infile)
        handle_error("Error: Input tree file doesn't find!");
    else{
        cout << "Successful open tree file!" << endl;
        //可以不用空格作为分隔符，一次读入一行
        while(getline(infile,line)) {
            strip(line);
            for(auto &c:line){
                if(c == '(')
                    store_flag = 1;
                if(store_flag && c == ';')
                    store_flag = 0;
                if(store_flag)
                    _tree_text += c;
            }
        }
        _tree_text += ';';
        infile.close();
    }
}

inline int Parse_Tree::_statistic_bracket_num(string &tree_text){
    int lef = 0, rig = 0;
    for(auto every_word = tree_text.begin(); every_word < tree_text.end(); every_word++){
        if(*every_word == '(') lef += 1;
        if(*every_word == ')') rig += 1;
    }
    return lef-rig;
}

inline string Parse_Tree::_make_tree_name(string &tree_text){
    int lef_ex_right_bracket_num = _statistic_bracket_num(tree_text);
    string tree_name;

    // 这个地方是利用右边一定是没有括号的，所以用右边把左边的括号去掉，只留下纯文本信息
    // 因为右边的括号可能藏在文本中，所以左边有可能还是以括号开头的，叶子节点是一定成立的
    tree_text = tree_text.substr(lef_ex_right_bracket_num, tree_text.size()-lef_ex_right_bracket_num);
    while(1){
        tree_name = random_string(3);
        if(!name_text_map.count(tree_name)){
            name_text_map[tree_name] = tree_text;
            break;
        }
    }
    return tree_name;
}

// 返回两个哈希，存储树上各个节点之间的父子关系
inline string Parse_Tree::_make_relationship(string &tree_text, string &father_tree_name){
    string child_tree_name = _make_tree_name(tree_text);
    child_parent_map[child_tree_name] = father_tree_name;
    parent_child_map[father_tree_name].push_back(child_tree_name);
    return child_tree_name;
}

// 根据树名称,检查树文本是否含有逗号确定是否是叶子节点
inline bool Parse_Tree::check_leaf_node(string &tree_name){
    string tree_text = name_text_map[tree_name];
    for(string::iterator i = tree_text.begin(); i < tree_text.end(); i++)
        if (*i == ',')
            return false;
    return true;
}


inline double Parse_Tree::give_tree_evolution_distance_father(string &tree_name){
    string tree_text = name_text_map[tree_name];
    string evol_dis = split(tree_text, ":").back();
    return stod(evol_dis);
}


// 这个函数也需要考虑树的模式，有区别的只有 >=3 有根树，但是 WHOLE 是根的单节点，因此没必要特殊处理
// <3 的有根树肯定是没问题
// >=3 的无根树 和 <3 的有根树 树的文本形式其实是一样的，不需要特殊处理
// 检验值和进化距离不一样，并不会给我们的程序带来误差波动
inline double Parse_Tree::give_test_result_self(string &tree_name){
    string tree_text = name_text_map[tree_name];
    smatch match_result;
    if(t_bootstrap && t_sh){
        // 特殊考虑叶子节点的情况
        if (check_leaf_node(tree_name))
            return 100;
        else{
            if(regex_search(tree_text, match_result, _test_double_pattern))
                return (stod(match_result[1]) + stod(match_result[2]))/2;
            else
                return 100;
        }
    }
    else{
        // 特殊考虑叶子节点的情况
        if(check_leaf_node(tree_name))
            return 100;
        else{
            if(regex_search(tree_text, match_result, _test_single_pattern))
                return stod(match_result[1]);
            else
                return 100;
        }
    }
}

// 这个函数不用考虑树的根，虽然 WHOLE 的时间不一定存在
// 需要查找 WHOLE 时间的时候，一定不是水平基因转移
inline double Parse_Tree::search_iqtree2_time(string &tree_name){
    string tree_text;
    smatch temp;
    tree_text = name_text_map[tree_name];
    regex_search(tree_text, temp, _time_pattern);
    return stod(temp[1]);
}

inline void Parse_Tree::judge_tree_root(){
    string text = split(name_text_map["WHOLE"], ":").back();
    int len = text.size();
    string new_text = text.substr(0,len);
    if(parent_child_map["WHOLE"].size() >= 3)
        // 根据冒号对树进行分割，如果最后一部分出现了括号，说明是三个(或者三个以上)节点无根树的形式
        if(new_text.find(')') != string::npos)
            cout << "Warning: the tree file you provide is an unroot tree, which may lead to uncertainty on our conclusion." << endl;
        else
            cout << "We detect the tree you provide is a root tree." << endl;
    else
        cout << "We detect the tree you provider is a root tree." << endl;
}


inline void Parse_Tree::print_child_tree(){
    cout << "CHILD TREE:" << endl;
    cout << "*******************************************************************************************" << endl;
    for(auto &it : parent_child_map){
        cout << name_text_map[it.first] << endl;
        for(auto & itchild : it.second)
            cout << name_text_map[itchild] << endl;
        cout << "*******************************************************************************************" << endl;
    }
}

inline void  Parse_Tree::print_tree_leaf(){
    cout << "LEAF NODE:" << endl;
    cout << "-------------------------------------------------------------------------------------------" << endl;
    for(auto &it : tree_leaf_map){
        cout << name_text_map[it.first] << endl;
        for(auto & itleaf : it.second)
            cout << itleaf << endl;
        cout << "-------------------------------------------------------------------------------------------" << endl;
    }
}
#endif

