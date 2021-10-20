#include "cloud.h"

// cloud 
cloud::cloud(vector<point>point_list, int data_num, int dimension, Comparator* comparator){
    this->point_list = point_list;
    this->data_num = data_num;
    this->dimension = dimension;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data_num);
    // cout<<"this is cloud!"<<endl;
}

void cloud::calculate_ef_xN(int N_index, vector<vector<int>>&e, vector<vector<int>>&f){
    /**
     * func:计算 ei = xi - ai, fi = Ni - bi
     * para: i -- the index of kd tree
     * *checked
     */ 
    for(int i = 0; i < data_num; ++i){
        // 理论上来说，f只要一列就好了，但是这里为了保持一致性，就和e规模一样，但是每一行都一样啦
        for(int j = 0; j < dimension; ++j){
            e[i][j] = point_list[i].data[j] - beaver_list[i][0];    // eij = xij - ai
            f[i][j] = kd_tree[N_index][i] - beaver_list[i][1];      // fij = Ni - bi
        }
    }
}
void cloud::calculate_ef_xN_square(vector<vector<int>>xN, vector<vector<int>>&e, vector<vector<int>>&f){
    // calculate (xi * Ni)^2 middle result e and f
    // eij = xNij - ai,     fij = xNij - bi
    // e, f has already been initialized
    for(int i = 0; i < data_num; ++i){
        for(int j = 0; j < dimension; ++j){
            e[i][j] = xN[i][j] - beaver_list[i][0];     // e = x - a
            f[i][j] = xN[i][j] - beaver_list[i][1];     // f = x - b
        }
    }
}

void cloud::calculate_ef_vari(vector<int>sum_xN, vector<vector<int>>&ef){
    // calculate (Σxi*Ni)^2 middle result e and f
    // first dimension is e, second dimension is f
    for(int i = 0; i < ef.size(); ++i){
        ef[i][0] = sum_xN[i] - beaver_list[i][0];   // e
        ef[i][1] = sum_xN[i] - beaver_list[i][1];   // f
    }
}

// cloud_one
cloud_one::cloud_one(vector<point>point_list, int data_num, int dimension, Comparator* comparator):
            cloud(point_list, data_num, dimension, comparator){
                // cout<<"this is cloud_one"<<endl;
            }

vector<vector<int>> cloud_one::calculate_xi_Ni(vector<vector<int>>&e, vector<vector<int>>&f){
    // xi*Ni_1 = fi * a1 + ei * b1 + c1
    vector<vector<int>>result(data_num, vector<int>(dimension,0));
    for(int i = 0; i < data_num; ++i){
        for(int j = 0; j < dimension; ++j){
            result[i][j] = f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<vector<int>> cloud_one::calculate_xi_Ni_square(vector<vector<int>>&e, vector<vector<int>>&f){
    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>>result(e);
    for(int i = 0; i < data_num; ++i){
        for(int j = 0; j < dimension; ++j){
            //? 看起来像是和上一个函数一毛一样（雀实，除了名字
            //? 应该是可以合并吧，但是蒜了，凑点代码显示工作量是不是
            result[i][j] = f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<int> cloud_one::calculate_sec_part(vector<vector<int>>ef){
    // ef是一个二维数组，rol1: e    rol2: f
    vector<int>result(ef.size(), 0);
    for(int i = 0; i < result.size(); ++i){
        result[i] = ef[i][1] * beaver_list[i][0] + ef[i][0] * beaver_list[i][1] + beaver_list[i][2];
    }
    return result;
}


// cloud_two 
cloud_two::cloud_two(vector<point>point_list, int data_num, int dimension, Comparator* comparator, vector<vector<int>>sorted_index):
            cloud(point_list, data_num, dimension, comparator){
                this->sorted_index = sorted_index;
                // cout<<"this is cloud two"<<endl;
            }

vector<vector<int>> cloud_two::calculate_xi_Ni(vector<vector<int>>&e, vector<vector<int>>&f){
    // cloud two calculate xi*Ni_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>>result(data_num, vector<int>(dimension,0));

    for(int i = 0; i < data_num; ++i)
        for(int j = 0; j < dimension; ++j)
            result[i][j] = e[i][j] * f[i][j] + f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];

    return result;        
}

vector<vector<int>> cloud_two::calculate_xi_Ni_square(vector<vector<int>>&e, vector<vector<int>>&f){
    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>>result(e);
    for(int i = 0; i < data_num; ++i){
        for(int j = 0; j < dimension; ++j){
            //? 看起来像是和上一个函数一毛一样（雀实，除了名字
            //? 应该是可以合并吧，但是蒜了，凑点代码显示工作量是不是
            result[i][j] = e[i][j] * f[i][j] + f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<int> cloud_two::calculate_sec_part(vector<vector<int>>ef){
    vector<int>result(ef.size(), 0);
    for(int i = 0; i < ef.size(); ++i)
        result[i] = ef[i][0] * ef[i][1] + ef[i][1] * beaver_list[i][0] + ef[i][0] * beaver_list[i][1] + beaver_list[i][2];

    return result;
}

