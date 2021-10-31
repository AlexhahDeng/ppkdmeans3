#pragma once

#include <iostream>
#include <ctime>
#include <vector>
#include <random>
#include <map>
#include <cmath>
#include <fstream>
#include <sstream>
#include "cloud.h"

using namespace std;

//! tools for myself

// 对每个维度的数据进行快速排序
void quick_sort(vector<point>& s, int dimension, int l, int r);

// 随机生成数据，但是不考虑能否正确聚类
void random(vector<point>& point_list, int n, int m);

// 读取需要聚类的数据
void read_data();

// 拆分数据
void divide_data(vector<point>point_list, int data_range, vector<point>&c1, vector<point>&c2);

// 对数据排序
vector<vector<int>> sort_data(vector<point>point_list);

// 生成beaver三元组
void generate_beaver_set(int n, int data_range, vector<vector<int>>&c1_list, vector<vector<int>>&c2_list);

// 固定化测试数据
void tmp_data(vector<point> &point_list);

// func: 用户初始化簇中心，确保点不重复--checked
void ini_clu_cen(cloud_one &c1, cloud_two &c2);

// func: 用户初始化候选中心和初始簇中心给root--checked
void ini_candidate_k(cloud_one &c1, cloud_two &c2);

// func: 预计算簇中心连乘的结果
void mul_clu_point_num(cloud_one& c1, cloud_two& c2);

