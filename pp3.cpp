#include <iostream>
// #include "cloud.h"
#include "func.h"
// #include "tools.h"

using namespace std;
/**
 * func： 构造kd tree
 */
void generate_kd_tree(cloud_one& c1, cloud_two& c2){
    // 从root开始，首先secret share N{1,1,1,...}
    vector<int>N(c1.data_num, 1);
    int num_point = c1.data_num;    //FIXME： 这是当前tree node中的点数
    c2.add_new_node(N, c1.kd_tree, c2.kd_tree);

    // 计算xi*Ni
    vector<vector<int>>e1(c1.data_num, vector<int>(c1.dimension, 0)), e2(e1);
    vector<vector<int>>f1(e1), f2(e1);

    // 计算xi*Ni的中间结果 e1,f1 and e2, f2
    c1.calculate_ef_xN(0, e1, f1);
    c2.calculate_ef_xN(0, e2, f2);

    //* 计算合并结果 e=e1+e2, f=f1+f2,那么就用e1和f1暂存中间结果e，f吧！
    for(int i = 0; i < c1.data_num; ++i){
        for(int j = 0; j < c1.dimension; ++j){
            e1[i][j] += e2[i][j];
            f1[i][j] += f2[i][j];
        } 
    }
    
    cout<<"计算xi*Ni..."<<endl;
    vector<vector<int>>xN1 = c1.calculate_xi_Ni(e1, f1);    //* checked
    vector<vector<int>>xN2 = c2.calculate_xi_Ni(e1, f1);

    // 接下来还得计算(xi*Ni)^2啊,mad
    // 1. 首先算中间结果e1,f1 and e2,f2
    c1.calculate_ef_xN_square(xN1, e1, f1);
    c2.calculate_ef_xN_square(xN2, e2, f2);

    // 2. 合并e，f，在e1,f1中暂存结果
    for(int i = 0; i < c1.data_num; ++i){
        for(int j = 0; j < c1.dimension; ++j){
            e1[i][j] += e2[i][j];
            f1[i][j] += f2[i][j];
        }
    }

    // 3. c1，c2分别计算(xi*Ni)^2
    cout<<"计算(xi*Ni)**2..."<<endl;
    vector<vector<int>>xN1_sq = c1.calculate_xi_Ni_square(e1, f1);  //* checked
    vector<vector<int>>xN2_sq = c2.calculate_xi_Ni_square(e1, f1);

    //* 计算方差了，uu们
    // 1. 首先计算(Σxi*Ni)^2， nΣ(xi*Ni)^2
    vector<int>sum_xN1(c1.dimension,0), sum_xN2(sum_xN1), sum_xN1_sq(sum_xN1), sum_xN2_sq(sum_xN1);

    for(int i = 0; i < c1.dimension; ++i){
        for(int j = 0; j < c1.data_num; ++j){
            sum_xN1[i] += xN1[j][i];
            sum_xN2[i] += xN2[j][i];
            sum_xN1_sq[i] += xN1_sq[j][i];
            sum_xN2_sq[i] += xN2_sq[j][i];
        }
        // nΣ(xi*Ni)^2 get， 改了改了
        // sum_xN1_sq[i] *= c1.data_num;
        // sum_xN2_sq[i] *= c2.data_num;
    }

    // 2. 再用一个beaver 三元组，算(Σxi*Ni)^2
    vector<vector<int>>ef1(c1.dimension, vector<int>(2)), ef2(ef1); // u1s1，我刚刚怎么没想到用一个二维数组存放ef，反正都在同一边！
    c1.calculate_ef_vari(sum_xN1, ef1);
    c2.calculate_ef_vari(sum_xN2, ef2);

    for(int i = 0; i < c1.dimension; ++i){  // ef1 = ef1 + ef2
        ef1[i][0] += ef2[i][0];
        ef1[i][1] += ef2[i][1];
    }

    vector<int>sec_part1 = c1.calculate_sec_part(ef1, num_point);
    vector<int>sec_part2 = c2.calculate_sec_part(ef1, num_point);

    // 3. 最后一步，组合起来！
    // 先在本地做减法，再加密组合成方差值
    for(int i = 0; i < c1.dimension; i++){
        sum_xN1_sq[i] -= sec_part1[i];
        sum_xN2_sq[i] -= sec_part2[i];
    }// 结果暂存在sum_xNi_sq中

    /* above has been checked */

    //TODO 要不还是留个接口出来，以后再处理把！
    // 1. c1和c2分别加密s1和s2
    vector<Ctxt>zero_one = c1.comparator->encrypt_variance(vector<int>{0,4369}, false); //FIXME：不可以直接加密1
    vector<Ctxt>enc_s1 = c1.comparator->encrypt_variance(sum_xN1_sq, true);
    vector<Ctxt>enc_s2 = c2.comparator->encrypt_variance(sum_xN2_sq, true);
    for(int i = 0; i < enc_s2.size(); ++i)  // 暂存在enc_s1中
        enc_s1[i] += enc_s2[i];


    // 2. c1 对密文{s1...sd}进行比较，求出最大值index，发给c2
    // Ctxt enc_index = c1.max_variance(enc_s1); //-->俺好怕，加起来以后超过模数范围咋整，再说把！

    // 3. c2 解密index，取 N 和 对应index的排序结果，划分，划分完了以后，秘密共享N_l N_r
    // c2.divide_data_set(enc_index, c1.kd_tree[0], 0, c1.data_num, c1); //* 介里第二个参数可以变化的啦

    return; 
}

void tmp_data(vector<point>&point_list){
    point_list.push_back({0, vector<int>{123, 43}});
    point_list.push_back({1, vector<int>{76, 987}});
    point_list.push_back({2, vector<int>{234, 76}});
    point_list.push_back({3, vector<int>{67, 567}});
    point_list.push_back({4, vector<int>{59, 473}});
}

int main(){
    Comparator* comparator = generate_comparator(true);
    int data_num = 5, dimension = 2;
    vector<point>point_list, c1_data, c2_data;
    tmp_data(point_list);
    // random(point_list, data_num, dimension);
    vector<vector<int>>sorted_index = sort_data(point_list);
    divide_data(point_list, 10000, c1_data, c2_data);

    cloud_one c1(c1_data, data_num, dimension, comparator);
    cloud_two c2(c2_data, data_num, dimension, comparator, sorted_index);

    generate_beaver_set(data_num, 100, c1.beaver_list, c2.beaver_list);

    // generate kd tree
    generate_kd_tree(c1, c2);
    // comparator->test_compare();
    return 0;
}