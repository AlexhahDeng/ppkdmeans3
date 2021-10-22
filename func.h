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
          