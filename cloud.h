#pragma once
#include <iostream>
#include <vector>
#include "tools.h"
using namespace std;

class cloud{
    public:
    int data_num;           // number of data record
    int dimension;          // dimension of a single point
    Comparator* comparator; // comparator for later comparison
    vector<vector<int>>beaver_list; // beaver triple list for multiplication
    vector<vector<int>>kd_tree;     // store data division info

    cloud(int data_num, int dimension, Comparator* comparator);

    vector<vector<int>> calculate_e_f(vector<int>N);
};

class cloud_one: public cloud{
    public:

};

class cloud_two: public cloud{

};