#pragma once
#include <iostream>
#include <vector>
#include "tools.h"

using namespace std;

struct point {
    int index;
    vector<int>data;
};

class cloud{
    public:
    int data_num;                   // number of data record
    int dimension;                  // dimension of a single point
    Comparator* comparator;         // comparator for later comparison
    vector<point>point_list;        // store data for secret sharing 
    vector<vector<int>>beaver_list; // beaver triple list for multiplication
    vector<vector<int>>kd_tree;     // store data division info

    cloud(vector<point>point_list, int data_num, int dimension, Comparator* comparator);
    
    // 计算xi * Ni的中间结果e,f
    void calculate_ef_xN(int N_index, vector<vector<int>>&e, vector<vector<int>>&f);

    // 计算(xi*Ni)^2 的中间结果e，f
    void calculate_ef_xN_square(vector<vector<int>>xN, vector<vector<int>>&e, vector<vector<int>>&f);

    // 计算(Σxi*Ni)^2的中间结果ef（存放在一个二维数组中
    void calculate_ef_vari(vector<int>sum_xN, vector<vector<int>>&ef);

};

class cloud_one: public cloud{
    public:
    //! 按道理来说有个shuffling rules，这里先不写（问题可能有点大，但我不想写

    cloud_one(vector<point>point_list, int data_num, int dimension, Comparator* comparator);

    // cloud one calculate xi*Ni_1 = f *a1 + e * b1 + c1
    vector<vector<int>> calculate_xi_Ni(vector<vector<int>>&e, vector<vector<int>>&f);

    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>> calculate_xi_Ni_square(vector<vector<int>>&e, vector<vector<int>>&f);

    // cloud_one calculate (Σxi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<int> calculate_sec_part(vector<vector<int>>ef);
};

class cloud_two: public cloud{
    public:
    vector<vector<int>>sorted_index;

    cloud_two(vector<point>point_list, int data_num, int dimension, Comparator* comparator, vector<vector<int>>sorted_index);

    // cloud two calculate xi*Ni_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>> calculate_xi_Ni(vector<vector<int>>&e, vector<vector<int>>&f);

    // cloud_two calculate (xi*Ni)^2_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>> calculate_xi_Ni_square(vector<vector<int>>&e, vector<vector<int>>&f);

    // cloud_two calculate (Σxi*Ni)^2_2 = e*f + f*a2 + e*b2 +c2
    vector<int> calculate_sec_part(vector<vector<int>>ef);
};