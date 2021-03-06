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
    long int node_point_num;          // number of points in current node
    vector<Ctxt> ctxt_node_point_num; // 虽然只有一个数，但是这样比较方便初始化
    vector<int> N;                    // secret share of division info, size = data_num
    vector<int> node_sum_x;           // Σxi*Ni, size = dimension
    vector<int> ptxt_node_min;        // 记录明文node中的最大最小值-->麻了，不ss了
    vector<int> ptxt_node_max;
    vector<Ctxt> ctxt_node_sum; // 记录密文和(非秘密共享值)
    vector<Ctxt> node_min;      // store min value of every dimension in curr node, size = dimension
    vector<Ctxt> node_max;      // store max value of every dimension in curr node, size = dimension
    vector<Ctxt> ctxt_candi_k;
    vector<int> candidate_k; // possible cluster belonging to

    bool isClustered = false; // 判断是否已经属于某个簇，true--则不进行任何操作，将child的属性设置为true；false则继续计算
};

class cloud
{
public:
    int k;         // number of clusters
    int data_num;  // number of data record
    int dimension; // dimension of a single point
    // Ctxt ctxt_one
    Comparator *comparator;          // comparator for later comparison
    vector<point> point_list;        // store data for secret sharing
    vector<vector<int>> beaver_list; // beaver triple list for multiplication, size = data_num×3
    vector<kd_node> kd_tree;         // store data division info

    vector<vector<int>> clu_cen; // size = k × dimension, no need to initialize
    vector<int> clu_point_num;   // len = k, update during iteration
    vector<int> mul_point_num;   // multiply the number of points in clusters

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
    vector<vector<int>> calculate_po_num_ef(vector<int> &v);

    // 计算距离的中间结果para的e，f
    vector<vector<long int>> calculate_dist_para_ef(int node_index, int k_index);

    // 计算距离结果的ef
    vector<vector<long int>> calculate_dist_res_ef(vector<long int> para);

    // 计算avg*Ni --> size= data_num*dimension
    void calculate_avg_N(vector<vector<int>> &e, vector<vector<int>> &f, vector<int> avg, int node_index);

    // 计算v*|z|的中间结果
    vector<vector<long int>> cal_vznum_ef(vector<long int> v, int k_index);

    // 计算自己的平方！
    vector<vector<long int>> cal_vec_square(vector<int> arr);

    vector<vector<long int>> cal_vec_square(vector<long int> arr);

    // 更新簇中心
    void update_clu_cen_ef(int node_index, vector<int> ptxt_mark, vector<vector<int>> &e, vector<vector<int>> &f);
};

class cloud_one : public cloud
{
public:
    vector<vector<Ctxt>> ctxt_clu_cen; // 密文的簇中心数据，方便在非叶节点部分，计算距离，进行prune操作

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
    vector<int> calculate_mul_final(vector<vector<int>> &ef);

    // cloud one calculate secret sharing distance parameters
    vector<long int> calculate_dist_para(vector<vector<long int>> ef);

    // cloud one calculate final distance result
    long int calculate_dist_res(vector<vector<long int>> ef);

    // cal v*|z| final res
    vector<long int> cal_vznum_final(vector<vector<long int>> ef);

    //cal square final
    vector<long int> cal_square_final(vector<vector<long int>> ef);

    // cal update clu cen final result
    vector<vector<int>> update_clu_cen_final(vector<vector<int>> &e, vector<vector<int>> &f);
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
    /**
     * @brief 添加新的节点
     * 
     * @param N 需要秘密拆分的位置信息数组，可以为负数，相加为01
     * @param point_num tree node包含的节点个数
     * @param c1_kdtree c1 维护的kd tree信息
     * @param c2_kdtree c2 维护的kd tree信息
     */
    void add_new_node(vector<int> N, int point_num, vector<kd_node> &c1_kdtree, vector<kd_node> &c2_kdtree);

    // 解密最大参数的下标
    int decrypt_index(Ctxt enc_var_index);

    // cloud two using ef to calculate final result e*f + f*α2 + e*b2 + c2
    vector<int> calculate_mul_final(vector<vector<int>> &ef);

    // cloud two calculate secret sharing distance parameters
    vector<long int> calculate_dist_para(vector<vector<long int>> ef);

    // cloud two calculate final distance result
    long int calculate_dist_res(vector<vector<long int>> ef);

    //
    vector<long int> cal_vznum_final(vector<vector<long int>> ef);

    vector<long int> cal_square_final(vector<vector<long int>> ef);

    // cal update clu cen final result
    vector<vector<int>> update_clu_cen_final(vector<vector<int>> &e, vector<vector<int>> &f);
};
