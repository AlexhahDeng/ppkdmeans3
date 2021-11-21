#pragma once

#include <iostream>
#include <ctime>
#include <vector>
#include <random>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include <time.h>
#include "cloud.h"

using namespace std;

//! tools for myself

// 对每个维度的数据进行快速排序
void quick_sort(vector<point>& s, int dimension, int l, int r);

// 随机生成数据，但是不考虑能否正确聚类
void random(vector<point>& point_list, int n, int m);

// 读取需要聚类的数据
vector<point> read_data(int data_num, int dimension);
// 拆分数据
void divide_data(vector<point>point_list, int data_range, vector<point>&c1, vector<point>&c2);

// 对数据排序
vector<vector<int>> sort_data(vector<point>point_list);

/**
 * @brief 生成beaver三元组
 * 
 * @param n beaver 三元组的个数 
 * @param data_range beaver 生成随机数模的范围
 * @param c1_list 赋给c1
 * @param c2_list 赋给c2
 */
void generate_beaver_set(int n, int data_range, vector<vector<int>>&c1_list, vector<vector<int>>&c2_list);

// 固定化测试数据
void tmp_data(vector<point> &point_list);

// func: 
/**
 * @brief 用户初始化簇中心，确保点不重复--checked
 * 
 * @param c1 存放了密文的初始簇中心以及秘密共享值
 * @param c2 仅存放秘密共享值
 */
void ini_clu_cen(cloud_one &c1, cloud_two &c2);

// func: 
/**
 * @brief 用户初始化候选中心和初始簇中心标识给root--checked
 * 
 * @param c1 初始每个簇都是候选{1,1,...}
 * @param c2 初始每个簇仅包含一个点-->不影响计算
 */
void ini_candidate_k(cloud_one &c1, cloud_two &c2);

/**
 * @brief 预计算簇中心连乘的结果，最后get{α1...αk,α}存入c1 and c2 *checked
 */
void mul_clu_point_num(cloud_one& c1, cloud_two& c2);

/**
 * @brief 向量内互不相同元素按顺序两两相乘，返回数组的size=k
 * 
 * @param c1 
 * @param c2 
 * @param node_index 该节点到簇的距离
 * @return vector<int> 到每个簇的距离, 除去不在候选中心的距离-->设置为0
 */
vector<int> cal_dist(cloud_one &c1, cloud_two &c2, int node_index, int& tot_can_num);
