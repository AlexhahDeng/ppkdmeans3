#include "cloud.h"

// cloud
cloud::cloud(vector<point> point_list, int data_num, int dimension, Comparator *comparator, int k)
{
    this->point_list = point_list;
    this->data_num = data_num;
    this->dimension = dimension;
    this->comparator = comparator;
    this->beaver_list = vector<vector<int>>(data_num);
    this->k = k;
    // this->clu_point_num = vector<int>(k);
    // this->clu_cen = vector<vector<int
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
    comparator->test_compare(1);

    return vec;
}

// v is single dimension
vector<vector<int>> cloud::calculate_po_num_ef(vector<int> &v)
{
    vector<vector<int>> result(v.size() / 2, vector<int>(2)); // size = v.size()/2 × 2
    for (int i = 0; i < v.size() / 2; ++i)
    {
        result[i][0] = v[2 * i] - beaver_list[i][0];     // e = x - a, x = α[2i]
        result[i][1] = v[2 * i + 1] - beaver_list[i][1]; // f = y - b, y = α[2i+1]
    }                                                    // 如果长度为偶数，最后一个数，不处理哦！
    return result;
}

vector<vector<int>> cloud::calculate_dist_res_ef(vector<int>para){
    // 0--Σc^2, 1--Σk^2, 2--Σci*ki, 3--α^2, 4--αj^2, 5--ααj，para的长度为6
    vector<vector<int>>result(3, vector<int>(2)); // 总共三项

    int i = 0;
    result[i][0] = para[0] - beaver_list[i][0]; // c^2 - a = e
    result[i][1] = para[3] - beaver_list[i][1]; // α^2 - b = f

    i++;
    result[i][0] = para[2] - beaver_list[i][0]; // Σci*ki - a = e
    result[i][1] = para[5] - beaver_list[i][1]; // ααj - b = f

    i++;
    result[i][0] = para[4] - beaver_list[i][0]; // αj^2 - a = e
    result[i][1] = para[1] - beaver_list[i][1]; // Σk^2 - b = f

    return result;
}


vector<vector<int>> cloud::calculate_dist_para_ef(int node_index, int k_index)
{
    vector<vector<int>> ef(dimension * 3 + 3, vector<int>(2));
    // 按照node中心{c0...cd},{k0...kd},{α, αj}的顺序依次排开
    int i = 0;
    while (i < dimension)
    {
        ef[i][0] = kd_tree[node_index].node_sum_x[i] / kd_tree[node_index].node_point_num - beaver_list[i][0]; // ci/n - a
        ef[i][1] = kd_tree[node_index].node_sum_x[i] / kd_tree[node_index].node_point_num - beaver_list[i][1]; // ci/n - b
        i++;
    } // 计算Ci×Ci

    while (i < dimension * 2)
    {
        ef[i][0] = clu_cen[k_index][i - dimension] - beaver_list[i][0]; // ki - a
        ef[i][1] = clu_cen[k_index][i - dimension] - beaver_list[i][1]; // ki - b
        ++i;
    } // 计算Ki×Ki

    while (i < dimension * 3)
    {
        ef[i][0] = kd_tree[node_index].node_sum_x[i - 2 * dimension] / kd_tree[node_index].node_point_num - beaver_list[i][0];  // ci - a
        ef[i][1] = clu_cen[k_index][i - 2 * dimension] - beaver_list[i][1];                                                     // ki - b    
        ++i;
    } // 计算 ci*ki

    ef[i][0] = mul_point_num[k] - beaver_list[i][0]; // α - a
    ef[i][1] = mul_point_num[k] - beaver_list[i][1]; // α - b
    ++i;

    ef[i][0] = mul_point_num[k_index] - beaver_list[i][0]; // αj - a
    ef[i][1] = mul_point_num[k_index] - beaver_list[i][1]; // αj - b
    i++;

    ef[i][0] = mul_point_num[k] - beaver_list[i][0];       // α - a
    ef[i][1] = mul_point_num[k_index] - beaver_list[i][1]; // αj - b

    return ef;
}

// cloud provides clu_cen, while node provide candidate_k
void cloud::exclude_clusters(int node_index){

}


void cloud::calculate_avg_N(vector<vector<int>>&e, vector<vector<int>>&f, vector<int>avg, int node_index)
{
    for(int i=0;i<data_num;i++){
        for(int j=0;j<dimension;j++){
            e[i][j] = kd_tree[node_index].N[i] - beaver_list[i][0]; // N[i]-a
            f[i][j] = avg[j] - beaver_list[i][1];                   // avg[j]-b
        }
    }
}

// cloud_one
cloud_one::cloud_one(vector<point> point_list, int data_num, int dimension, Comparator *comparator, int k) : cloud(point_list, data_num, dimension, comparator, k)
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

vector<int> cloud_one::calculate_mul_final(vector<vector<int>> &ef)
{
    vector<int> result(ef.size());
    for (int i = 0; i < ef.size(); ++i)
    {
        result[i] = beaver_list[i][0] * ef[i][1] + beaver_list[i][1] * ef[i][0] + beaver_list[i][2];
    }
    return result;
}

vector<int> cloud_one::calculate_dist_para(vector<vector<int>> ef)
{
    vector<int> para(ef.size());
    for (int i = 0; i < ef.size(); i++)
    {
        para[i] = ef[i][0] * beaver_list[i][1] + ef[i][1] * beaver_list[i][0] + beaver_list[i][2];
    }
   vector<int>result(6);
    int sum_c_square = 0, sum_k_square = 0, c_mul_k = 0, i = 0;
    // Σc^2
    while(i < dimension)
        sum_c_square += para[i++];
    result[0] = sum_c_square;

    // Σk^2
    while(i < 2*dimension)
        sum_k_square += para[i++];
    result[1] = sum_k_square;

    // Σc*k
    while(i < 3*dimension)
        c_mul_k += para[i++];
    result[2] = c_mul_k;

    result[3] = para[i++]; // α^2
    result[4] = para[i++]; // αj^2
    result[5] = para[i++]; // α*αj

    return result;
}

int cloud_one::calculate_dist_res(vector<vector<int>>ef){
    vector<int>para(ef.size());
    for(int i = 0; i < ef.size(); i++){
        para[i] = ef[i][0]*beaver_list[i][1] + ef[i][1]*beaver_list[i][0] + beaver_list[i][2];
    }
    return para[0] - 2*para[1] + para[2];
}


// cloud_two
cloud_two::cloud_two(vector<point> point_list, int data_num, int dimension, Comparator *comparator, vector<vector<int>> sorted_index, int k) : cloud(point_list, data_num, dimension, comparator, k)
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
    int max_var_index = comparator->decrypt_index(enc_index); // 解密下标值
    vector<int> n_min(dimension), n_max(n_min);                    // 记录每个维度的最大最小值
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
                    n_min[j] = c1.point_list[index].data[j] + this->point_list[index].data[j];
                }
                if (count_mm[j] == curr_data_num)
                {
                    n_max[j] = c1.point_list[index].data[j] + this->point_list[index].data[j];
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
    c1.kd_tree[N_index].node_min = comparator->encrypt_vector(n_min);
    c1.kd_tree[N_index].node_max = comparator->encrypt_vector(n_max);


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
    vector<Ctxt>ctxt_node_point=comparator->encrypt_vector(vector<int>{point_num});

    c1_kdtree.push_back({point_num, ctxt_node_point, N1, vector<int>(dimension), vector<Ctxt>(), vector<Ctxt>(), vector<Ctxt>(), vector<Ctxt>(), vector<int>(k)});
    c2_kdtree.push_back({point_num, ctxt_node_point, N2, vector<int>(dimension), vector<Ctxt>(), vector<Ctxt>(), vector<Ctxt>(), vector<Ctxt>(), vector<int>(k)});
}

vector<int> cloud_two::calculate_mul_final(vector<vector<int>> &ef)
{
    vector<int> result(ef.size());
    for (int i = 0; i < ef.size(); ++i)
    {
        result[i] = ef[i][0] * ef[i][1] + beaver_list[i][0] * ef[i][1] + beaver_list[i][1] * ef[i][0] + beaver_list[i][2];
    }
    return result;
}

vector<int> cloud_two::calculate_dist_para(vector<vector<int>> ef)
{
    vector<int> para(ef.size());
    for (int i = 0; i < ef.size(); i++)
    {
        para[i] = ef[i][1]*ef[i][0] + ef[i][0] * beaver_list[i][1] + ef[i][1] * beaver_list[i][0] + beaver_list[i][2];
    }

    vector<int>result(6);
    int sum_c_square = 0, sum_k_square = 0, c_mul_k = 0, i = 0;
    // Σc^2
    while(i < dimension)
        sum_c_square += para[i++];
    result[0] = sum_c_square;

    // Σk^2
    while(i < 2*dimension)
        sum_k_square += para[i++];
    result[1] = sum_k_square;

    // Σc*k
    while(i < 3*dimension)
        c_mul_k += para[i++];
    result[2] = c_mul_k;

    result[3] = para[i++]; // α^2
    result[4] = para[i++]; // αj^2
    result[5] = para[i++]; // α*αj

    return result;
}


int cloud_two::calculate_dist_res(vector<vector<int>>ef){
    vector<int>para(ef.size());
    for(int i = 0; i < ef.size(); i++){
        para[i] = ef[i][0]*ef[i][1] + ef[i][0]*beaver_list[i][1] + ef[i][1]*beaver_list[i][0] + beaver_list[i][2];
    }
    return para[0] - 2*para[1] + para[2];
}