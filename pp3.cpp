#include <iostream>
#include "cloud.h"
#include "func.h"
#include "tools.h"

using namespace std;
/**
 * func： 构造kd tree
 */
void generate_kd_tree(cloud_one& c1, cloud_two& c2){
    // 从root开始，首先secret share N{1,1,1,...}
    vector<int>N(c1.data_num, 1);
    secret_share_N(N, c1, c2);

    // 计算xi*Ni
    vector<vector<int>>e1(c1.data_num, vector<int>(c1.dimension, 0)), e2(e1);
    vector<vector<int>>f1(e1), f2(e1);

    c1.calculate_ef_xN(0, e1, f1);
    c2.calculate_ef_xN(0, e2, f2);

    //* 计算合并结果，e,f,那么就用e1和f1暂存中间结果e，f吧！
    for(int i = 0; i < c1.data_num; ++i){
        for(int j = 0; j < c1.dimension; ++j){
            e1[i][j] += e2[i][j];
            f1[i][j] += f2[i][j];
        } 
    }

    vector<vector<int>>xN1 = c1.calculate_xi_Ni(e1, f1);
    vector<vector<int>>xN2 = c2.calculate_xi_Ni(e1, f1);

    // 接下来还得计算(xi*Ni)^2啊mad
    // 1. 首先算中间结果
    c1.calculate_ef_xN_square(xN1, e1, f1);
    c2.calculate_ef_xN_square(xN2, e2, f2);

    // 2. 合并e，f，在e1中暂存结果
    for(int i = 0; i < c1.data_num; ++i){
        for(int j = 0; j < c1.dimension; ++j){
            e1[i][j] += e2[i][j];
            f1[i][j] += f2[i][j];
        }
    }

    // 3. c1，c2分别计算(xi*Ni)^2
    vector<vector<int>>xN1_sq = c1.calculate_xi_Ni_square(e1, f1);
    vector<vector<int>>xN2_sq = c2.calculate_xi_Ni_square(e1, f1);
    
    return; 
}
int main(){
    Comparator* comparator = generate_comparator(true);

    // // initialization 
    // int data_num = 5, dimension = 2;
    // vector<point>point_list, c1_data, c2_data;
    // random(point_list, data_num, dimension);
    // vector<vector<int>>sorted_index = sort_data(point_list);
    // divide_data(point_list, 10000, c1_data, c2_data);

    // cloud_one c1(c1_data, data_num, dimension, comparator);
    // cloud_two c2(c2_data, data_num, dimension, comparator, sorted_index);

    // generate_beaver_set(data_num, 100, c1.beaver_list, c2.beaver_list);

    // // generate kd tree
    // generate_kd_tree(c1, c2)      ;
    comparator->sort();
    return 0;
}