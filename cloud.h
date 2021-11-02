#pragma once
#include <iostream>
#include <vector>
#include "tools.h"
#include <ctime>
#include <random>
using namespace std;

struct point
{
    int index;
    vector<int> data;
};

struct kd_node
{
    int node_point_num;               // number of points in current node
    vector<int> N;                    // secret share of division info, size = data_num
    vector<int> node_sum_x;           // Σxi*Ni, size = dimension
    vector<vector<int>> node_min_max; // store max and min value of every dimension in curr node, size = dimension * 2
    vector<int> candidate_k;          // possible cluster belonging to
};

class cloud
{
public:
    int k;                           // number of clusters
    int data_num;                    // number of data record
    int dimension;                   // dimension of a single point
    Comparator *comparator;          // comparator for later comparison
    vector<point> point_list;        // store data for secret sharing
    vector<vector<int>> beaver_list; // beaver triple list for multiplication, size = data_num×3
    vector<kd_node> kd_tree;         // store data division info

    vector<vector<int>> clu_cen;     // size = k × dimension, no need to initialize
    vector<int>clu_point_num;        // len = k, update during iteration
    vector<int>mul_point_num;        // multiply the number of points in clusters

    cloud(vector<point> point_list, int data_num, int dimension, Comparator *comparator, int k);

    // 计算xi * Ni的中间结果e,f
    void calculate_ef_xN(int N_index, vector<vector<int>> &e, vector<vector<int>> &f);

    // 计算(xi*Ni)^2 的中间结果e，f
    void calculate_ef_xN_square(vector<vector<int>> xN, vector<vector<int>> &e, vector<vector<int>> &f);

    // 计算(Σxi*Ni)^2的中间结果ef（存放在一个二维数组中
    void calculate_ef_vari(vector<int> sum_xN, vector<vector<int>> &ef);

    // 加密输入的数据
    vector<Ctxt> encrypt_variance(vector<int> vari);

    // 计算(αi*αi+1)的中间结果e，f
    vector<vector<int>> calculate_po_num_ef(vector<int>&v);

    // 计算距离的中间结果e，f
    vector<vector<int>> calculate_dist_ef(int node_index, int k_index);
};

class cloud_one : public cloud
{
public:
    //! 按道理来说有个shuffling rules，这里先不写（问题可能有点大，但我不想写

    cloud_one(vector<point> point_list, int data_num, int dimension, Comparator *comparator, int k);

    // cloud one calculate xi*Ni_1 = f *a1 + e * b1 + c1
    vector<vector<int>> calculate_xi_Ni(vector<vector<int>> &e, vector<vector<int>> &f);

    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>> calculate_xi_Ni_square(vector<vector<int>> &e, vector<vector<int>> &f);

    // cloud_one calculate (Σxi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<int> calculate_sec_part(vector<vector<int>> ef, int n);

    // cloud_one figure out index of maximum variance
    Ctxt max_variance(vector<Ctxt> enc_variance, vector<Ctxt> zero_one);

    // cloud one using ef to get final result f*a1 + e*b1 + c1
    vector<int> calculate_mul_final(vector<vector<int>>& ef);

    // cloud one calculate secret sharing distance parameters
    vector<int> calculate_dist_para(vector<int>ef);

};

class cloud_two : public cloud
{
public:
    vector<vector<int>> sorted_index; // size = dimension×data_num

    cloud_two(vector<point> point_list, int data_num, int dimension, Comparator *comparator, vector<vector<int>> sorted_index, int k);

    // cloud two calculate xi*Ni_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>> calculate_xi_Ni(vector<vector<int>> &e, vector<vector<int>> &f);

    // cloud_two calculate (xi*Ni)^2_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>> calculate_xi_Ni_square(vector<vector<int>> &e, vector<vector<int>> &f);

    // cloud_two calculate (Σxi*Ni)^2_2 = e*f + f*a2 + e*b2 +c2
    vector<int> calculate_sec_part(vector<vector<int>> ef, int n);

    /**
     * cloud two decrypt enc_cipher, divide data set using N, index and sorted result
     * enc_index: encrypted index for dimension
     * vector<int>N1: shuffled secret share N, used to recover N
     * N_index: index of N used for cloud two to get respect N2
     * tot_p: number of points in this tree node
     */
    void divide_data_set(Ctxt enc_index, cloud_one &c1, vector<int> &N, int curr_data_num, int N_index);

    // 生成秘密共享位置信息
    void add_new_node(vector<int> N, int point_num, vector<kd_node> &c1_kdtree, vector<kd_node> &c2_kdtree);

    // 解密最大参数的下标
    int decrypt_index(Ctxt enc_var_index);

    // cloud two using ef to calculate final result e*f + f*α2 + e*b2 + c2
    vector<int> calculate_mul_final(vector<vector<int>>& ef);

    // cloud one calculate secret sharing distance
    vector<int> calculate_dist_res(vector<int>ef);

    // cloud two calculate secret sharing distance parameters
    vector<int> calculate_dist_para(vector<int> ef);


};
