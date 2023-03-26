// Created by zhangshengbin on 2021/7/20.
// 献给毕生所爱
#ifndef _OTHER_F_H
#define _OTHER_F_H

#include <string>
#include <vector>
#include <map>

#include <cctype>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
#include <fstream>

using namespace std;

// 用来给各个子树随机生成名字的函数
inline string random_string(size_t);

// 去除头尾空格、换行和制表符
inline string strip(string &);

// 根据分隔符对字符串进行分割，返回向量
inline vector<string> split(string const &, string del = " ");

// 修改变量名，把不是数字、字母和下划线的部分改为数字、字母和下划线
inline void change_variable_name(string &name);
// 把修改的变量名改回去
inline string change_name_back(string name);

// 利用mmap把一个文件读到内存里，快速获取文件行数
inline char* map_file(const char* fname, size_t& length);

// 给出一个整数是几位数
inline int give_num_figure(int num);

// 给出各个聚类的名字，根据文件按行来命名的
inline string give_OG_name(int OG_num, int total_figure);

// 寻找最久远子树的进化距离
inline double find_min_map_value(map<string, double> &my_map);
// 寻找最久远的子树的名称
inline void find_min_map_key(map<string, double> &my_map, vector<string> &min_key, vector<string> &temp_tree_vec);

// kd树取中位数用到的函数，向 0 取整
inline double get_median(vector<double> vec);

// 读取蛋白汇总文件中的蛋白名称和序列关系
inline void read_summary_protein_file(string &protein_file_name, map<string, string> &prot_name_sequence_map);

// 调换各个派发线程哈希的前后关系，使任务数量一致
inline void move_map_order(map<string, vector<string> > &ori_map, map<string, vector<string> > &target_map, int threads, int grainsize);

// 输出错误信息
inline void handle_error(const char* msg);
static map<string, string> name_change_map;

inline string random_string(size_t length)
{
    auto randchar = []() -> char
    {
        const char charset[] =
                "0123456789"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[rand() % max_index];
    };
    string str(length,0);
    generate_n(str.begin(), length, randchar);
    return str;
}

inline vector<string> split(string const &s, string del){
    int start = 0;
    int end = s.find(del);
    string token;
    vector<string> res;
    while (end != -1){
        token = s.substr(start, end - start);
        res.push_back(token);
        start = end + del.size();
        end = s.find(del, start);
    }
    token = s.substr(start, end - start);
    res.push_back(token);
    return res;
}

inline string strip(string &inpt){
    auto start_it = inpt.begin();
    auto end_it = inpt.rbegin();
    while (isspace(*start_it))
        ++start_it;
    while (isspace(*end_it))
        ++end_it;
    return string(start_it, end_it.base());
}

inline void change_variable_name(string &name){
    string temp(name);
    int flag = 0;
    for(auto &c:name)
        if (!isalnum(c) && c != '_'){
            c = '_';
            flag = 1;
        }
    if(flag == 1)
        name_change_map[temp] = name;
}

inline string change_name_back(string name){
    for(auto &it : name_change_map)
        if(it.second == name)
            return it.first;
    return name;
}

inline void handle_error(const char* msg){
    perror(msg);
    exit(255);
}

inline char* map_file(const char* fname, size_t& length)
{
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("Error: Input orthogroups file doesn't find!");
    // obtain file size
    struct stat sb{};
    if (fstat(fd, &sb) == -1)
        handle_error("fstat error");
    length = sb.st_size;
    char* addr = static_cast<char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap error");
    // 在合适的地方使用 munmap 清除内存
    return addr;
}

inline int give_num_figure(int num){
    int figure = 1;
    while(num/10 >= 1){
        figure++;
        num = num/10;
    }
    return figure;
}

inline string give_OG_name(int OG_num, int total_figure){
    int OG_figure = give_num_figure(OG_num);
    int dif = total_figure - OG_figure, zero_n = 0;
    string zero = "", final_name;
    while(zero_n < dif){
        zero += "0";
        zero_n++;
    }
    final_name = "OG" + zero + to_string(OG_num);
    return  final_name;
}

inline double find_min_map_value(map<string, double> &my_map){
    double temp_min = 10000000;
    for(auto &eve_k:my_map){
        auto temp_val = eve_k.second;
        if(temp_min > temp_val)
            temp_min = temp_val;
    }
    return temp_min;
}

inline void find_min_map_key(map<string, double> &my_map, vector<string> &min_key, vector<string> &temp_tree_vec){
    double min = find_min_map_value(my_map);
    for(auto eve_k = my_map.begin(); eve_k != my_map.end(); eve_k++) {
        // 浮点数比大小是真行，直接用等号是不行得有误差，这里设置一个误差线
        if ((eve_k->second) - min < 0.0001)
            min_key.push_back(eve_k->first);
        else
            temp_tree_vec.push_back(eve_k->first);
    }
}

// 返回一个 vector 中的中位数
inline double get_median(vector<double> vec){
    if(vec.size() != 0){
        sort(vec.begin(), vec.end());
        int half = floor(vec.size()/2);
        return vec[half];
    }
    else
        return 0;
}

// 读取所有蛋白的蛋白文件，并生成一个字典准备进行输出
inline void read_summary_protein_file(string &protein_file_name, map<string, string> &prot_name_sequence_map){
    ifstream infile(protein_file_name);
    string line;
    int first_protein_flag = 0;

    if(!infile)
        handle_error("Error: Input protein file doesn't find!");

    else{
        cout << "Successful open protein sequence file!" << endl;
        string every_prot_name, protein_sequence;
        while(getline(infile, line)){

            if(line[0] == '>'){
                // getline 自动丢弃换行符
                if(first_protein_flag == 0){
                    first_protein_flag = 1;
                    protein_sequence = "";
                }

                else{
                    prot_name_sequence_map[every_prot_name] = protein_sequence;
                    protein_sequence = "";
                }

                // 这个是去掉＞号的
                every_prot_name = split(line, " ")[0];
                int len = every_prot_name.size();
                every_prot_name = every_prot_name.substr(1,len-1);
                change_variable_name(every_prot_name);

            }
            else{
                protein_sequence += line;
                protein_sequence += '\n';
            }
        }
        prot_name_sequence_map[every_prot_name] = protein_sequence;
    }
}

inline void move_map_order(map<string, vector<string> > &ori_map, map<string, vector<string> > &target_map, int threads, int grainsize){
    auto iter_begin = ori_map.begin();
    auto iter_temp = iter_begin;
    for(int k = 0; k < grainsize; k++)
    {
        iter_begin = iter_temp;
        for(int i = 0; i < threads-1; i++){
            target_map[iter_begin->first] = iter_begin->second;
            if(i != threads-2)
                for(int j = 0; j < grainsize; j++)
                    iter_begin ++;
        }
        iter_temp++;
    }
    iter_begin ++;
    for(auto iter = iter_begin; iter != ori_map.end(); iter++)
        target_map[iter->first] = iter->second;
}


#endif