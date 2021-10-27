#include "cloud.h"

// cloud
cloud::cloud(vector<point> point_list, int data_num, int dimension, Comparator *comparator)
{
    this->point_list = point_list;
    this->data_num = data_num;
    this->dimension = dimension;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data_num);
    // this->kd_tree = vector<vector<kd_node>>();
    // cout<<"this is cloud!"<<endl;
}

void cloud::calculate_ef_xN(int N_index, vector<vector<int>> &e, vector<vector<int>> &f)
{
    /**
     * func:计算 ei = xi - ai, fi = Ni - bi
     * para: i -- the index of kd tree
     * *checked
     */
    for (int i = 0; i < data_num; ++i)
    {
        // 理论上来说，f只要一列就好了，但是这里为了保持一致性，就和e规模一样，但是每一行都一样啦
        for (int j = 0; j < dimension; ++j)
        {
            e[i][j] = point_list[i].data[j] - beaver_list[i][0]; // eij = xij - ai
            f[i][j] = kd_tree[N_index].N[i] - beaver_list[i][1]; // fij = Ni - bi
        }
    }
}

void cloud::calculate_ef_xN_square(vector<vector<int>> xN, vector<vector<int>> &e, vector<vector<int>> &f)
{
    // calculate (xi * Ni)^2 middle result e and f
    // eij = xNij - ai,     fij = xNij - bi
    // e, f has already been initialized
    for (int i = 0; i < data_num; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            e[i][j] = xN[i][j] - beaver_list[i][0]; // e = x - a
            f[i][j] = xN[i][j] - beaver_list[i][1]; // f = x - b
        }
    }
}

void cloud::calculate_ef_vari(vector<int> sum_xN, vector<vector<int>> &ef)
{
    // calculate (Σxi*Ni)^2 middle result e and f
    // first dimension is e, second dimension is f
    for (int i = 0; i < ef.size(); ++i)
    {
        ef[i][0] = sum_xN[i] - beaver_list[i][0]; // e
        ef[i][1] = sum_xN[i] - beaver_list[i][1]; // f
    }
}

vector<Ctxt> cloud::encrypt_variance(vector<int> vari)
{
    /**
     * func: 明文方差-->密文，争取能够逐数字加密，即，one number-->one cipher
     * ? 实在不行的话，那就multiple same number --> one cipher
     * ? or 如果会shift的话，可以用那个多轮比较的办法，但是不好说
     * ? why not lib min function? 因为非常固定化，并且结果是方差，不是01的结果
     * TODO 冲啊！
     */

    vector<Ctxt> vec;
    comparator->test_compare();

    return vec;
}

// cloud_one
cloud_one::cloud_one(vector<point> point_list, int data_num, int dimension, Comparator *comparator) : cloud(point_list, data_num, dimension, comparator)
{
    // cout<<"this is cloud_one"<<endl;
}

vector<vector<int>> cloud_one::calculate_xi_Ni(vector<vector<int>> &e, vector<vector<int>> &f)
{
    // xi*Ni_1 = fi * a1 + ei * b1 + c1
    vector<vector<int>> result(data_num, vector<int>(dimension, 0));
    for (int i = 0; i < data_num; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            result[i][j] = f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<vector<int>> cloud_one::calculate_xi_Ni_square(vector<vector<int>> &e, vector<vector<int>> &f)
{
    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>> result(e);
    for (int i = 0; i < data_num; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            //? 看起来像是和上一个函数一毛一样（雀实，除了名字
            //? 应该是可以合并吧，但是蒜了，凑点代码显示工作量是不是
            result[i][j] = f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<int> cloud_one::calculate_sec_part(vector<vector<int>> ef, int n)
{
    // ef是一个二维数组，rol1: e    rol2: f
    vector<int> result(ef.size(), 0);
    for (int i = 0; i < result.size(); ++i)
    {
        result[i] = int((ef[i][1] * beaver_list[i][0] + ef[i][0] * beaver_list[i][1] + beaver_list[i][2]) / n);
    }
    return result;
}

Ctxt cloud_one::max_variance(vector<Ctxt> enc_variance, vector<Ctxt> zero_one)
{
    /**
     * func: 比较密文方差的大小，得到最大值的密文index
     * * 参考comparator部分的代码，自己写哦，没法直接调用库函数，因为俺把加密和比较的过程拆开了
     * ! 要特别注意乘法深度，和噪声的问题，可以手动实现一下 “两两比较”
     */

    return zero_one[0];
}

// cloud_two
cloud_two::cloud_two(vector<point> point_list, int data_num, int dimension, Comparator *comparator, vector<vector<int>> sorted_index) : cloud(point_list, data_num, dimension, comparator)
{
    this->sorted_index = sorted_index;
    // cout<<"this is cloud two"<<endl;
}

vector<vector<int>> cloud_two::calculate_xi_Ni(vector<vector<int>> &e, vector<vector<int>> &f)
{
    // cloud two calculate xi*Ni_2 = e*f + f*a2 + e*b2 + c2
    vector<vector<int>> result(data_num, vector<int>(dimension, 0));

    for (int i = 0; i < data_num; ++i)
        for (int j = 0; j < dimension; ++j)
            result[i][j] = e[i][j] * f[i][j] + f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];

    return result;
}

vector<vector<int>> cloud_two::calculate_xi_Ni_square(vector<vector<int>> &e, vector<vector<int>> &f)
{
    // cloud_two calculate (xi*Ni)^2_1 = f*a1 + e*b1 + c1
    vector<vector<int>> result(e);
    for (int i = 0; i < data_num; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            //? 看起来像是和上一个函数一毛一样（雀实，除了名字
            //? 应该是可以合并吧，但是蒜了，凑点代码显示工作量是不是
            result[i][j] = e[i][j] * f[i][j] + f[i][j] * beaver_list[i][0] + e[i][j] * beaver_list[i][1] + beaver_list[i][2];
        }
    }
    return result;
}

vector<int> cloud_two::calculate_sec_part(vector<vector<int>> ef, int n)
{
    // n is the number of points in this tree node
    vector<int> result(ef.size(), 0);
    for (int i = 0; i < ef.size(); ++i)
        result[i] = int((ef[i][0] * ef[i][1] + ef[i][1] * beaver_list[i][0] + ef[i][0] * beaver_list[i][1] + beaver_list[i][2]) / n);

    return result;
}

int cloud_two::decrypt_index(Ctxt enc_var_index)
{

    return 1;
}

void cloud_two::divide_data_set(Ctxt enc_index, cloud_one &c1, vector<int> &N, int curr_data_num, int N_index)
{
    /**
     * func: 1. merge N1 and N2 as N(noted that although N is plaintext, its order has been shuffled)
     *       2. decrypt index and obtain according sorted index
     *       3. divide result Nl and Nr, send Nl1, Nr1 to cloud 1
     * points:
     *       1. how to decrypt and decode data?
     *       2. 因为主函数也要算N，干脆直接在那里算完！
     */

    // TODO decrypt and decode ciphertext index
    int max_var_index = 0; // get_index(enc_index);   //-->解密的函数有现成的，但是怎么decode需要自己写
    vector<vector<int>> node_min_max1(dimension, vector<int>(2)), node_min_max2(node_min_max1);
    int count = 0;
    vector<int> count_mm(dimension, 0);
    vector<int> Nl(N), Nr(N);
    srand(time(NULL));

    for (int i = 0; i < data_num; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            int index = sorted_index[j][i];
            if (N[index])
            {
                count_mm[j]++;
                if (count_mm[j] == 1)
                {
                    node_min_max1[j][0] = rand() % dimension;
                    node_min_max2[j][0] = index - node_min_max1[j][0];
                }
                if (count_mm[j] == curr_data_num)
                {
                    node_min_max1[j][1] = rand() % dimension;
                    node_min_max2[j][1] = index - node_min_max1[j][1];
                }
            }
        } // checked统计node中每个维度数据的最小、最大值
        if (N[sorted_index[max_var_index][i]])
        {
            if (count < curr_data_num / 2)
                Nr[sorted_index[max_var_index][i]] = 0;
            else
                Nl[sorted_index[max_var_index][i]] = 0;
            count++;
        }
    }
    // add min_max to node info
    c1.kd_tree[N_index].node_min_max = node_min_max1;
    this->kd_tree[N_index].node_min_max = node_min_max2;

    // add new node to kd_tree
    add_new_node(Nl, curr_data_num / 2, c1.kd_tree, this->kd_tree);
    add_new_node(Nr, curr_data_num - curr_data_num / 2, c1.kd_tree, this->kd_tree);
}

void cloud_two::add_new_node(vector<int> N, int point_num, vector<kd_node> &c1_kdtree, vector<kd_node> &c2_kdtree)
{
    /**
     * func: 将由01构成的N数组，划分为由01构成的N1，N2两个点集，并分别存入c1和c2中
     * *checked
     */
    vector<int> N1(N), N2(N);
    srand(time(NULL));
    for (int i = 0; i < N.size(); ++i)
    {
        N1[i] = rand() % 2;
        N2[i] = N[i] - N1[i];
        // 不是做的加减法，是做的异或！
        // UPDATE: 只能做加减法，不然后续乘法部分会出错哦！
    }

    c1_kdtree.push_back({point_num, N1, vector<int>(dimension), vector<vector<int>>(dimension, vector<int>(2))});
    c2_kdtree.push_back({point_num, N2, vector<int>(dimension), vector<vector<int>>(dimension, vector<int>(2))});
}
