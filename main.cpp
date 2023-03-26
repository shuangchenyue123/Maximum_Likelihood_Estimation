// 弄完这个之后用itol做另外两幅图吧，先把这一部分结果彻底统计好。
// 然后写python多线程的脚本，一步导出COG结果我们基本上就差不多完成了。

/*
 2022.2.21任务 目前还是不行，还需要统计一下每个阶段发生水平基因转移的 ortho 有多大概率是比我们提供的数字大的。
 我们当时的假设就是计算发生HGT ortho 的 p，不发生的没考虑，所以这里没问题
*/
#include <iostream>
#include <cmath>
#include <unordered_map>
#include "ortho_class.h"

void parse_events_file(string& file_name, map<int, unordered_map<string, int>> &ortho_HGT_present_dic){
    string line;
    ifstream infile(file_name.c_str());
    while(getline(infile, line)){
        if(line.substr(0,4) == "gene")
            continue;
        else{
            string temp = strip(line);
            vector<string> temp_vec = split(temp,"\t");
            if(temp_vec.back() != "cluster1")
                continue;
            double HGT_time = stod(temp_vec[1]);
            string ortho_name = temp_vec[3];
            if(HGT_time < 100){
                if(!ortho_HGT_present_dic[0].count(ortho_name))
                    ortho_HGT_present_dic[0].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[0][ortho_name]++;
            }

            else if(HGT_time >= 100 and HGT_time < 200){
                if(!ortho_HGT_present_dic[100].count(ortho_name))
                    ortho_HGT_present_dic[100].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[100][ortho_name]++;
            }

            else if(HGT_time >= 200 and HGT_time < 300){
                if(!ortho_HGT_present_dic[200].count(ortho_name))
                    ortho_HGT_present_dic[200].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[200][ortho_name]++;
            }

            else if(HGT_time >= 300 and HGT_time < 400){
                if(!ortho_HGT_present_dic[300].count(ortho_name))
                    ortho_HGT_present_dic[300].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[300][ortho_name]++;
            }

            else if(HGT_time >= 400 and HGT_time < 500){
                if(!ortho_HGT_present_dic[400].count(ortho_name))
                    ortho_HGT_present_dic[400].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[400][ortho_name]++;
            }

            else if(HGT_time >= 500 and HGT_time < 600){
                if(!ortho_HGT_present_dic[500].count(ortho_name))
                    ortho_HGT_present_dic[500].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[500][ortho_name]++;
            }

            else if(HGT_time >= 600 and HGT_time < 700){
                if(!ortho_HGT_present_dic[600].count(ortho_name))
                    ortho_HGT_present_dic[600].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[600][ortho_name]++;
            }

            else if(HGT_time >= 700 and HGT_time < 800){
                if(!ortho_HGT_present_dic[700].count(ortho_name))
                    ortho_HGT_present_dic[700].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[700][ortho_name]++;
            }

            else if(HGT_time >= 800 and HGT_time < 900){
                if(!ortho_HGT_present_dic[800].count(ortho_name))
                    ortho_HGT_present_dic[800].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[800][ortho_name]++;
            }

            else if(HGT_time >= 900 and HGT_time < 1000){
                if(!ortho_HGT_present_dic[900].count(ortho_name))
                    ortho_HGT_present_dic[900].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[900][ortho_name]++;
            }

            else if(HGT_time >= 1000 and HGT_time < 1100){
                if(!ortho_HGT_present_dic[1000].count(ortho_name))
                    ortho_HGT_present_dic[1000].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1000][ortho_name]++;
            }

            else if(HGT_time >= 1100 and HGT_time < 1200){
                if(!ortho_HGT_present_dic[1100].count(ortho_name))
                    ortho_HGT_present_dic[1100].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1100][ortho_name]++;
            }

            else if(HGT_time >= 1200 and HGT_time < 1300){
                if(!ortho_HGT_present_dic[1200].count(ortho_name))
                    ortho_HGT_present_dic[1200].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1200][ortho_name]++;
            }

            else if(HGT_time >= 1300 and HGT_time < 1400){
                if(!ortho_HGT_present_dic[1300].count(ortho_name))
                    ortho_HGT_present_dic[1300].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1300][ortho_name]++;
            }

            else if(HGT_time >= 1400 and HGT_time < 1500){
                if(!ortho_HGT_present_dic[1400].count(ortho_name))
                    ortho_HGT_present_dic[1400].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1400][ortho_name]++;
            }

            else if(HGT_time >= 1500 and HGT_time < 1600){
                if(!ortho_HGT_present_dic[1500].count(ortho_name))
                    ortho_HGT_present_dic[1500].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1500][ortho_name]++;
            }

            else if(HGT_time >= 1600 and HGT_time < 1700){
                if(!ortho_HGT_present_dic[1600].count(ortho_name))
                    ortho_HGT_present_dic[1600].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1600][ortho_name]++;
            }

            else if(HGT_time >= 1700 and HGT_time < 1800){
                if(!ortho_HGT_present_dic[1700].count(ortho_name))
                    ortho_HGT_present_dic[1700].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1700][ortho_name]++;
            }

            else if(HGT_time >= 1800 and HGT_time < 1900){
                if(!ortho_HGT_present_dic[1800].count(ortho_name))
                    ortho_HGT_present_dic[1800].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1800][ortho_name]++;
            }

            else if(HGT_time >= 1900 and HGT_time < 2000){
                if(!ortho_HGT_present_dic[1900].count(ortho_name))
                    ortho_HGT_present_dic[1900].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[1900][ortho_name]++;
            }

            else if(HGT_time >= 2000 and HGT_time < 2100){
                if(!ortho_HGT_present_dic[2000].count(ortho_name))
                    ortho_HGT_present_dic[2000].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2000][ortho_name]++;
            }

            else if(HGT_time >= 2100 and HGT_time < 2200){
                if(!ortho_HGT_present_dic[2100].count(ortho_name))
                    ortho_HGT_present_dic[2100].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2100][ortho_name]++;
            }

            else if(HGT_time >= 2200 and HGT_time < 2300){
                if(!ortho_HGT_present_dic[2200].count(ortho_name))
                    ortho_HGT_present_dic[2200].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2200][ortho_name]++;
            }

            else if(HGT_time >= 2300 and HGT_time < 2400){
                if(!ortho_HGT_present_dic[2300].count(ortho_name))
                    ortho_HGT_present_dic[2300].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2300][ortho_name]++;
            }

            else if(HGT_time >= 2400 and HGT_time < 2500){
                if(!ortho_HGT_present_dic[2400].count(ortho_name))
                    ortho_HGT_present_dic[2400].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2400][ortho_name]++;
            }

            else if(HGT_time >= 2500 and HGT_time < 2600){
                if(!ortho_HGT_present_dic[2500].count(ortho_name))
                    ortho_HGT_present_dic[2500].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2500][ortho_name]++;
            }

            else if(HGT_time >= 2600 and HGT_time < 2700){
                if(!ortho_HGT_present_dic[2600].count(ortho_name))
                    ortho_HGT_present_dic[2600].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2600][ortho_name]++;
            }

            else if(HGT_time >= 2700 and HGT_time < 2800){
                if(!ortho_HGT_present_dic[2700].count(ortho_name))
                    ortho_HGT_present_dic[2700].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2700][ortho_name]++;
            }

            else if(HGT_time >= 2800 and HGT_time < 2900){
                if(!ortho_HGT_present_dic[2800].count(ortho_name))
                    ortho_HGT_present_dic[2800].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2800][ortho_name]++;
            }

            else if(HGT_time >= 2900 and HGT_time < 3000){
                if(!ortho_HGT_present_dic[2900].count(ortho_name))
                    ortho_HGT_present_dic[2900].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[2900][ortho_name]++;
            }

            else if(HGT_time >= 3000 and HGT_time < 3100){
                if(!ortho_HGT_present_dic[3000].count(ortho_name))
                    ortho_HGT_present_dic[3000].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[3000][ortho_name]++;
            }


            else if(HGT_time >= 3100 and HGT_time < 3200){
                if(!ortho_HGT_present_dic[3100].count(ortho_name))
                    ortho_HGT_present_dic[3100].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[3100][ortho_name]++;
            }

            else if(HGT_time >= 3200 and HGT_time < 3300){
                if(!ortho_HGT_present_dic[3200].count(ortho_name))
                    ortho_HGT_present_dic[3200].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[3200][ortho_name]++;
            }


            else if(HGT_time >= 3300 and HGT_time < 3400){
                if(!ortho_HGT_present_dic[3300].count(ortho_name))
                    ortho_HGT_present_dic[3300].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[3300][ortho_name]++;
            }

            else if(HGT_time >= 3400 and HGT_time < 3500){
                if(!ortho_HGT_present_dic[3400].count(ortho_name))
                    ortho_HGT_present_dic[3400].insert(make_pair(ortho_name, 1));
                else
                    ortho_HGT_present_dic[3400][ortho_name]++;
            }

        }
    }
}

double get_P(double max_P, double min, double HGT_num_total, vector<double> &beita_vec){
    double first_sum = 0, P;
    while(1){
        double second_sum = 0;
        P = (max_P + min)/2;
        for(auto &j:beita_vec){
            second_sum += j/(1-j*P);
        }
        first_sum = HGT_num_total/P;
        if(abs(first_sum - second_sum) < 0.01)
                break;
        else {
            if (first_sum - second_sum > 0)
                min = P;
            else
                max_P = P;
        }
    }
    return P;
}

// 统计以下各个时间段的物种都有多少个
// 做法上应该是先把树解析了，然后直接遍历一遍，遍历的时候计算时间，统计各个时间段物种有多少
// 已知根节点在4000万万年之前了，这一点不需要纠结
int main(int argc,char *argv[]) {
    string tree_file_n = argv[2], tree_root = "WHOLE", orthofile = argv[3], eventsfile = argv[4];
    double min_t = stod(argv[1]);
    map<int, double> time_species_map;
    map<int, unordered_set<string>> time_node_map;
    map<int, double> time_ortho_map;

    // 下面这map用来解析events.tsv
    map<int, unordered_map<string, int>> ortho_HGT_present_dic;


    // 初始化字典，用前边界作为字典的键
    for(int i = 0; i < 3500; i+=100) {
        time_species_map.insert(make_pair(i, 0));
        time_ortho_map.insert(make_pair(1,0));
    }


    // 解析树文件
    Parse_Tree tree(tree_file_n, true, true, false);
    tree.read_tree();
    // 这个地方前面虽然是下划线开头，实际上是个public变量，当时是为了防止重名
    tree.parse_tree(tree._tree_text, tree_root);
    tree.judge_tree_root();
    tree.find_tree_leaf(tree_root);
    tree.make_tree_attribution_map(tree_root);
    tree.oat = tree.tree_evolution_map[tree.leaf_name_to_code_map[tree.tree_leaf_map[tree_root][0]]];

    // 把阈值算出来
    sliding_mode sl_matrix(tree, min_t);
    sl_matrix.calculate_every_node_threshold();




    // 好了，现在开始进行统计, tree_evolution_map里有各个节点的进化值
    // 里面没有根节点，不会搞出real_time为0，根节点也不需要加进去，因为自己和自己没啥好转的
    for(auto &i:tree.tree_evolution_map){
        double real_time = 4000 - i.second;
        string node_name = i.first;

        // 处理一下计算误差
        if(real_time < 0.01 and real_time > -0.01)
            real_time = 0;

        if(real_time < 100) {
            time_species_map[0]++;
            time_node_map[0].insert(node_name);
        }

        if(real_time >= 100 and real_time < 200) {
            time_species_map[100]++;
            time_node_map[100].insert(node_name);
        }

        else if(real_time >= 200 and real_time < 300) {
            time_species_map[200]++;
            time_node_map[200].insert(node_name);
        }

        else if(real_time >= 300 and real_time < 400) {
            time_species_map[300]++;
            time_node_map[300].insert(node_name);
        }

        else if(real_time >= 400 and real_time < 500) {
            time_species_map[400]++;
            time_node_map[400].insert(node_name);
        }

        else if(real_time >= 500 and real_time < 600) {
            time_species_map[500]++;
            time_node_map[500].insert(node_name);
        }

        else if(real_time >= 600 and real_time < 700) {
            time_species_map[600]++;
            time_node_map[600].insert(node_name);
        }

        else if(real_time >= 700 and real_time < 800) {
            time_species_map[700]++;
            time_node_map[700].insert(node_name);
        }

        else if(real_time >= 800 and real_time < 900) {
            time_species_map[800]++;
            time_node_map[800].insert(node_name);
        }

        else if(real_time >= 900 and real_time < 1000) {
            time_species_map[900]++;
            time_node_map[900].insert(node_name);
        }

        else if(real_time >= 1000 and real_time < 1100) {
            time_species_map[1000]++;
            time_node_map[1000].insert(node_name);
        }

        else if(real_time >= 1100 and real_time < 1200) {
            time_species_map[1100]++;
            time_node_map[1100].insert(node_name);
        }

        else if(real_time >= 1200 and real_time < 1300) {
            time_species_map[1200]++;
            time_node_map[1200].insert(node_name);
        }

        else if(real_time >= 1300 and real_time < 1400) {
            time_species_map[1300]++;
            time_node_map[1300].insert(node_name);
        }

        else if(real_time >= 1400 and real_time < 1500) {
            time_species_map[1400]++;
            time_node_map[1400].insert(node_name);
        }

        else if(real_time >= 1500 and real_time < 1600) {
            time_species_map[1500]++;
            time_node_map[1500].insert(node_name);
        }

        else if(real_time >= 1600 and real_time < 1700) {
            time_species_map[1600]++;
            time_node_map[1600].insert(node_name);
        }

        else if(real_time >= 1700 and real_time < 1800) {
            time_species_map[1700]++;
            time_node_map[1700].insert(node_name);
        }

        else if(real_time >= 1800 and real_time < 1900) {
            time_species_map[1800]++;
            time_node_map[1800].insert(node_name);
        }

        else if(real_time >= 1900 and real_time < 2000) {
            time_species_map[1900]++;
            time_node_map[1900].insert(node_name);
        }

        else if(real_time >= 2000 and real_time < 2100) {
            time_species_map[2000]++;
            time_node_map[2000].insert(node_name);
        }

        else if(real_time >= 2100 and real_time < 2200) {
            time_species_map[2100]++;
            time_node_map[2100].insert(node_name);
        }

        else if(real_time >= 2200 and real_time < 2300) {
            time_species_map[2200]++;
            time_node_map[2200].insert(node_name);
        }

        else if(real_time >= 2300 and real_time < 2400) {
            time_species_map[2300]++;
            time_node_map[2300].insert(node_name);
        }

        else if(real_time >= 2400 and real_time < 2500) {
            time_species_map[2400]++;
            time_node_map[2400].insert(node_name);
        }

        else if(real_time >= 2500 and real_time < 2600) {
            time_species_map[2500]++;
            time_node_map[2500].insert(node_name);
        }

        else if(real_time >= 2600 and real_time < 2700) {
            time_species_map[2600]++;
            time_node_map[2600].insert(node_name);
        }

        else if(real_time >= 2700 and real_time < 2800) {
            time_species_map[2700]++;
            time_node_map[2700].insert(node_name);
        }

        else if(real_time >= 2800 and real_time < 2900) {
            time_species_map[2800]++;
            time_node_map[2800].insert(node_name);
        }

        else if(real_time >= 2900 and real_time < 3000) {
            time_species_map[2900]++;
            time_node_map[2900].insert(node_name);
        }

        else if(real_time >= 3000 and real_time < 3100) {
            time_species_map[3000]++;
            time_node_map[3000].insert(node_name);
        }

        else if(real_time >= 3100 and real_time < 3200) {
            time_species_map[3100]++;
            time_node_map[3100].insert(node_name);
        }

        else if(real_time >= 3200 and real_time < 3300) {
            time_species_map[3200]++;
            time_node_map[3200].insert(node_name);
        }

        else if(real_time >= 3300 and real_time < 3400) {
            time_species_map[3300]++;
            time_node_map[3300].insert(node_name);
        }
        else if(real_time >= 3400 and real_time < 3500) {
            time_species_map[3400]++;
            time_node_map[3400].insert(node_name);
        }

        // 到3500可以了，因为HGT就到3500，在往上没有了

    }

    // 好了现在把orthofile给读取了
    OthoHGT ortho(orthofile ,tree, true, min_t, sl_matrix, time_node_map);
    ortho.read_file();
    ortho.statistic_HGT_genes_paraller(15);

    // 解析Events.file文件
    parse_events_file(eventsfile, ortho_HGT_present_dic);

    // 下面应该可以汇总解析了
    // 输出相应的信息到csv文件中
    ofstream out("Time_period_species.tsv");
    out << "time" << '\t' << "HGT_possibility" << endl;
    for(auto &time_sec:time_species_map) {
        out << time_sec.first << "~" << time_sec.first+100 << '\t';
        // time_sec_total_HGT就是x的和，用K来表示β的和，这两个都是一个ortho只有一个，α是多个系数算出来的，max_P记录了这段时间区间里P可以取的最大值, sample_num是k值
        double max_P = 1, HGT_num_total = 0, second_sum = 0, first_sum = 0;
        vector<double> beita_vec;
        // 处理时间区间内的各个ortho
        for(auto &every_ortho:ortho_HGT_present_dic[time_sec.first]){
            double ortho_node = 0, beita = 0;
            for(auto &every_node:time_node_map[time_sec.first])
                if(ortho.ortho_node_present[every_ortho.first].find(every_node) != ortho.ortho_node_present[every_ortho.first].end())
                    ortho_node++;

            // 如果发生了HGT，必然有2个或以上的node是含有这个ortho的，否则不可能
            int ortho_HGT_num = every_ortho.second;
            ortho_node = ortho_node - ortho_HGT_num;
            double P_contain = ortho_node/time_sec.second, P_not = 1 - P_contain, HGT_num;
            HGT_num_total += ortho_HGT_num;

            for(HGT_num = 1; HGT_num <= ortho_HGT_num; HGT_num++){
                max_P > 1/(time_sec.second*P_contain*P_not) ? max_P = 1/(time_sec.second*P_contain*P_not) : max_P = max_P;
                // 每一轮水平基因转移都更新一次
                P_contain = (ortho_node+1)/time_sec.second, P_not = 1 - P_contain;
            }
            max_P > 1/(time_sec.second*P_contain*P_not) ? max_P = 1/(time_sec.second*P_contain*P_not) : max_P = max_P;
            beita_vec.push_back(time_sec.second*P_contain*P_not);
            // k 代表每一个ortho得系数，HGT代表指数，K代表一个时间区段的总系数
        }
        max_P = max_P*0.999;
        for(auto &beita:beita_vec){
            second_sum += beita/(1-beita*max_P);
        }
        first_sum = HGT_num_total/max_P;
        if(first_sum - second_sum >= 0)
            out << max_P/ortho.time_ortho_map[time_sec.first].size() << endl;
        else{
            //求极值点
            double P = get_P(max_P, 0, HGT_num_total, beita_vec);
            out << P/ortho.time_ortho_map[time_sec.first].size() << endl;
        }
    }

    out.close();

}
