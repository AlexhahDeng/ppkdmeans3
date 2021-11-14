//
//                       _oo0oo_
//                      o8888888o
//                      88" . "88
//                      (| -_- |)
//                      0\  =  /0
//                    ___/`---'\___
//                  .' \\|     |// '.
//                 / \\|||  :  |||// \
//                / _||||| -:- |||||- \
//               |   | \\\  -  /// |   |
//               | \_|  ''\---/''  |_/ |
//               \  .-\__  '-'  ___/-. /
//             ___'. .'  /--.--\  `. .'___
//          ."" '<  `.___\_<|>_/___.' >' "".
//         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
//         \  \ `_.   \_ __\ /__ _/   .-` /  /
//     =====`-.____`.___ \_____/___.-`___.-'=====
//                       `=---='
//
//
//     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//               佛祖保佑         实验顺利
//
//
//
#include <iostream>
// #include "cloud.h"
#include "my_tools.h"
// #include "tools.h"

using namespace std;
// global para
// vector<Ctxt>zero_one = c1.comparator->encrypt_vector(vector<int>{0,4369}, false); //FIXME：要根据整数分解的位数灵活变化

/**
 * func： 构造kd tree
 */
void generate_kd_tree(cloud_one &c1, cloud_two &c2)
{
    // 从root开始，首先secret share N{1,1,1,...}
    vector<int> N(c1.data_num, 1);
    c2.add_new_node(N, c1.data_num, c1.kd_tree, c2.kd_tree);
    int N_index = 0;

    while (true)
    {
        // 计算当前tree node中包含点个数，组合N
        int num_point = 0;
        for (int i = 0; i < c1.data_num; ++i)
        {
            N[i] = c1.kd_tree[N_index].N[i] + c2.kd_tree[N_index].N[i]; // 还原秘密共享分解值
            if (N[i])
                num_point++;
        }
        if (num_point == 1)
            break;

        // 计算xi*Ni
        // 1. 初始化e1，e2 和 f1，f2，二维数组，data_num*dimension
        vector<vector<int>> e1(c1.data_num, vector<int>(c1.dimension, 0)), e2(e1);
        vector<vector<int>> f1(e1), f2(e1);

        // 2. 计算xi*Ni的中间结果 e1,f1 and e2, f2
        c1.calculate_ef_xN(N_index, e1, f1);
        c2.calculate_ef_xN(N_index, e2, f2);

        // 3. 计算合并结果 e=e1+e2, f=f1+f2,那么就用e1和f1暂存中间结果e，f吧！
        for (int i = 0; i < c1.data_num; ++i)
        {
            for (int j = 0; j < c1.dimension; ++j)
            {
                e1[i][j] += e2[i][j];
                f1[i][j] += f2[i][j];
            }
        }

        // 4. xN1 = f1*a1+e1*b1+c1  xN2=e2*f2+f2*a2+e2*b2+c2
        cout << "计算xi*Ni..." << endl;
        vector<vector<int>> xN1 = c1.calculate_xi_Ni(e1, f1); //* checked
        vector<vector<int>> xN2 = c2.calculate_xi_Ni(e1, f1);

        // 接下来还得计算(xi*Ni)^2啊,mad
        // 1. 首先算中间结果e1,f1 and e2,f2
        c1.calculate_ef_xN_square(xN1, e1, f1);
        c2.calculate_ef_xN_square(xN2, e2, f2);

        // 2. 合并e，f，在e1,f1中暂存结果
        for (int i = 0; i < c1.data_num; ++i)
        {
            for (int j = 0; j < c1.dimension; ++j)
            {
                e1[i][j] += e2[i][j];
                f1[i][j] += f2[i][j];
            }
        }

        // 3. c1，c2分别计算(xi*Ni)^2
        cout << "计算(xi*Ni)**2..." << endl;
        vector<vector<int>> xN1_sq = c1.calculate_xi_Ni_square(e1, f1); //* checked
        vector<vector<int>> xN2_sq = c2.calculate_xi_Ni_square(e1, f1);

        //* 计算方差了，uu们
        // 1. 首先计算(Σxi*Ni)^2， Σ(xi*Ni)^2
        vector<int> sum_xN1(c1.dimension, 0), sum_xN2(sum_xN1), sum_xN1_sq(sum_xN1), sum_xN2_sq(sum_xN1);

        for (int i = 0; i < c1.dimension; ++i)
        {
            for (int j = 0; j < c1.data_num; ++j)
            {
                sum_xN1[i] += xN1[j][i];
                sum_xN2[i] += xN2[j][i];
                sum_xN1_sq[i] += xN1_sq[j][i];
                sum_xN2_sq[i] += xN2_sq[j][i];
            }
            // nΣ(xi*Ni)^2 get， 改了改了
            // sum_xN1_sq[i] *= c1.data_num;
            // sum_xN2_sq[i] *= c2.data_num;
        }

        //! 将(Σxi*Ni)补充到kd_node中
        c1.kd_tree[N_index].node_sum_x = sum_xN1;
        c2.kd_tree[N_index].node_sum_x = sum_xN2;

        // 2. 再用一个beaver 三元组，算(Σxi*Ni)^2
        vector<vector<int>> ef1(c1.dimension, vector<int>(2)), ef2(ef1); // u1s1，我刚刚怎么没想到用一个二维数组存放ef，反正都在同一边！
        c1.calculate_ef_vari(sum_xN1, ef1);
        c2.calculate_ef_vari(sum_xN2, ef2);

        for (int i = 0; i < c1.dimension; ++i)
        { // ef1 = ef1 + ef2
            ef1[i][0] += ef2[i][0];
            ef1[i][1] += ef2[i][1];
        }

        vector<int> sec_part1 = c1.calculate_sec_part(ef1, num_point);
        vector<int> sec_part2 = c2.calculate_sec_part(ef1, num_point);

        // 3. 最后一步，组合起来！
        // 先在本地做减法，再加密组合成方差值
        vector<int> index(c1.dimension); // 因为后面获取最大方差，还需要index的密文
        for (int i = 0; i < c1.dimension; i++)
        {
            sum_xN1_sq[i] -= sec_part1[i];
            sum_xN2_sq[i] -= sec_part2[i];

            sum_xN1_sq[i] += sum_xN2_sq[i];
        } // 结果暂存在sum_xNi_sq中

        /* above has been checked */

        //! 目前先直接合并秘密共享的方差，然后加密-->就不用解决负数的问题了
        vector<Ctxt> enc_value = c1.comparator->encrypt_vector(sum_xN1_sq);

        // 2. c1 对密文{s1...sd}进行比较，求出最大值index，发给c2
        Ctxt ctxt_one = c1.comparator->gen_ctxt_one();
        Ctxt max_index = c1.comparator->max_variance(enc_value, ctxt_one); //! 没实现，就先默认第一个维度最大把
        c1.comparator->decrypt_index(max_index);

        // 3. c2 解密index，取 N 和 对应index的排序结果，划分，划分完了以后，秘密共享N_l N_r
        c2.divide_data_set(max_index, c1, N, num_point, N_index);
        N_index++;
        break;
    }

    return;
}

// func: 聚类过程
void filtering(cloud_one &c1, cloud_two &c2)
{

    ini_clu_cen(c1, c2);     // 初始化簇中心（同时加密存放在c1中），存放在cloud中
    ini_candidate_k(c1, c2); // 初始化候选集和簇中点数，存放在kd_node中

    while (true)
    { // 一轮迭代

        mul_clu_point_num(c1, c2);
        // 预计算簇中心连乘的结果，每一轮迭代只计算一次
        //* clu_point_num已经重置为0

        vector<vector<int>> new_clu_cen(c1.k, vector<int>(c1.dimension, 0));
        // 这里可以把密文结果累加，存储到一起，迭代完再进行ss划分成两份
        //? 存在的问题是-->中间累加的结果可能会超过加密的范围，所以我们简化问题
        // 直接存储明文结果把！

        vector<int> new_clu_point_num(c1.k);// 包含的点数同理，存储明文

        for (int i = 0; i < c1.kd_tree.size(); i++) // 遍历kd tree 所有节点
        {
            if (c1.kd_tree[i].isClustered){
                c1.kd_tree[2*i+1].isClustered = true;
                c1.kd_tree[2*i+2].isClustered = true;
                continue;
            } // 如果已经被划分了则不再继续划分,并且标识子节点不再划分

            vector<vector<int>> dist = cal_dist(c1, c2, i);
            // 计算中心到每个簇中心的距离, size = 2 x k -->其实也可以直接合并后加密啦……emmm也不太影响嘛
            
            int tot_candidate_k = 0;
            // 记录当前node包含多少个候选中心
            for (int j = 0; j < c1.k; j++) // 先直接合并，存到维度0中
            {
                int n = c1.kd_tree[i].candidate_k[j] + c2.kd_tree[i].candidate_k[j];
                // 合并n，候选为1，非候选为0

                tot_candidate_k += n;
                
                
                dist[0][j] = (dist[0][j] + dist[1][j]) * n;
                // 本来应该是要用ss乘法的，懒得实现了，麻了
            }

            vector<Ctxt> ctxt_dist = c1.comparator->encrypt_vector(dist[0]);
            // 加密距离-->可能面临超过范围的问题

            Ctxt ctxt_one = c1.comparator->gen_ctxt_one(); //目前没什么用，但是懒得改函数了
            vector<Ctxt> min_dist_index = c1.comparator->min_dist(ctxt_dist, ctxt_one);
            // 处理0影响的方法
            // 对每个距离密文增加一个近似最大值，0就变成了最大值，但是其他的会被模，相对大小不改变

            // 根据是否为叶子节点分不同情况处理
            if (c1.kd_tree[i].node_point_num > 2) // 非叶子节点
            {
                // 计算z*-->距离最近的簇
                // 构造密文0
                Ctxt ctxt_zero = ctxt_one;
                ctxt_zero *= 0l;
                Ctxt less_than = ctxt_one;
                vector<Ctxt> k_closest(c1.dimension, ctxt_zero);

                for (int j = 0; j < c1.dimension; j++)
                {
                    for (int k = 0; k < c1.k; k++)
                    {
                        Ctxt ctxt_value = min_dist_index[k];
                        ctxt_value *= long(c1.clu_cen[k][j] + c2.clu_cen[k][j]);
                        k_closest[j] += ctxt_value; //-->可能会超过范围哦！-->不会吧，本来就在范围内啊
                    }
                }

                // 遍历所有的簇中心，prune掉距离更远的
                vector<Ctxt> v(c1.dimension, ctxt_one);
                for (int j = 0; j < c1.k; j++)
                {
                    for (int k = 0; k < c1.dimension; k++)
                    {
                        Ctxt u = c1.ctxt_clu_cen[j][k];
                        u -= k_closest[k];                               // u[k] = z[k] - z*[k]
                        c1.comparator->compare(less_than, u, ctxt_zero); // u[i] < 0?

                        // 解密结果tmd
                        long int com_res = c1.comparator->dec_compare_res(less_than);
                        Ctxt mid_res = c1.kd_tree[i].node_min[k];
                        mid_res *= com_res; // LT(u[i],0)*(c_min)
                        v[k] = mid_res;

                        mid_res = c1.kd_tree[i].node_max[k];
                        mid_res *= (1 - com_res); //(1-LT(u[i],0))*(c_max)
                        v[k] += mid_res;
                    }// 首先计算当前簇对应的u

                    Ctxt ctxt_d1 = ctxt_zero, ctxt_d2 = ctxt_zero;
                    for (int k = 0; k < c1.dimension; k++)
                    {
                        Ctxt curr = k_closest[k];
                        curr -= v[k]; // (cloest_ki - vi)
                        curr *= curr; // (cloest_ki - vi)^2
                        ctxt_d1 += curr;

                        curr = c1.ctxt_clu_cen[j][k];
                        curr -= v[k]; // (ki - vi)
                        curr *= curr; // (ki - vi)**2
                        ctxt_d2 += curr;
                    }// 计算dist(v,z*), dist(v,z)

                    c1.comparator->compare(less_than, ctxt_d1, ctxt_d2); // dist(v,z*)<dist(v,z)?
                    long res = c1.comparator->dec_compare_res(less_than);
                    for (int k = 0; k < c1.dimension; k++)
                    {
                        c1.kd_tree[i].candidate_k[j] = res;
                        c2.kd_tree[i].candidate_k[j] = 0;

                        c1.kd_tree[2 * i + 1].candidate_k[j] = res;
                        c2.kd_tree[2 * i + 1].candidate_k[j] = 0; // 心中随机即可

                        c1.kd_tree[2 * i + 2].candidate_k[j] = res;
                        c2.kd_tree[2 * i + 2].candidate_k[j] = 0; // 心中随机即可
                    }// 修改当前node和子node的候选中心

                    tot_candidate_k -= res;
                    // 如果簇被prune，候选中心数目减少,当然论文里面咱们不能这么做
                }
                // 判断候选candidate是否只包含一个簇
                if(tot_candidate_k == 1){
                    //* 这里的步骤基本上和下面叶子节点的步骤一致
                    vector<long> node_cen(c1.dimension);
                    // 合并点集中心结果
                    for (int j = 0; j < c1.dimension; j++)
                        node_cen[j] = (c1.kd_tree[i].node_sum_x[j] + c2.kd_tree[i].node_sum_x[j]) / c1.kd_tree[i].node_point_num;
                    // 这里需要做除法，因为在算距离的地方做的是除法

                    vector<Ctxt> ctxt_node_cen = c1.comparator->encrypt_vector(node_cen);
                    // 加密点集中心

                    for (int j = 0; j < c1.k; j++)
                    {
                        for (int k = 0; k < c1.dimension; k++)
                        {
                            Ctxt curr = ctxt_node_cen[k];
                            curr *= min_dist_index[j];
                            // node每一个维度与对应k的min_dist结果相乘, be like {c1...cd}*0/1

                            int value = c2.comparator->decrypt_index(curr);
                            // 按理来说，这里应该先做分割再decrypt，蒜了，nevermind
                            new_clu_cen[j][k] += value;
                        }
                        // 首先将node中包含的点个数与min_dist结果相乘
                        Ctxt ctxt_point_num = min_dist_index[j];
                        ctxt_point_num *= c1.kd_tree[i].node_point_num;

                        // 秘密 共享 拆分-->那么直接解密把！
                        int res = c2.comparator->decrypt_index(ctxt_point_num);

                        new_clu_point_num[j] += res;
                    } // 将集合中的点拆分，分别存入c1，c2

                    // 标记当前和子节点信息
                    c1.kd_tree[i].isClustered = true;
                    c1.kd_tree[2*i+1].isClustered = true;
                    c1.kd_tree[2*i+2].isClustered = true;

                }// 更新簇信息，并且标识当前node和子node都不再进行聚类


            }
            else // 叶子节点
            {
                vector<long> node_cen(c1.dimension);
                // 合并点集中心结果
                for (int j = 0; j < c1.dimension; j++)
                    node_cen[j] = (c1.kd_tree[i].node_sum_x[j] + c2.kd_tree[i].node_sum_x[j]) / c1.kd_tree[i].node_point_num;
                // 这里需要做除法，因为在算距离的地方做的是除法

                vector<Ctxt> ctxt_node_cen = c1.comparator->encrypt_vector(node_cen);
                // 加密点集中心

                for (int j = 0; j < c1.k; j++)
                {
                    for (int k = 0; k < c1.dimension; k++)
                    {
                        Ctxt curr = ctxt_node_cen[k];
                        curr *= min_dist_index[j];
                        // node每一个维度与对应k的min_dist结果相乘, be like {c1...cd}*0/1

                        int value = c2.comparator->decrypt_index(curr);
                        // 按理来说，这里应该先做分割再decrypt，蒜了，nevermind
                        new_clu_cen[j][k] += value;
                    }
                    // 首先将node中包含的点个数与min_dist结果相乘
                    Ctxt ctxt_point_num = min_dist_index[j];
                    ctxt_point_num *= c1.kd_tree[i].node_point_num;

                    // 秘密 共享 拆分-->那么直接解密把！
                    int res = c2.comparator->decrypt_index(ctxt_point_num);

                    new_clu_point_num[j] += res;
                } // 将集合中的点拆分，分别存入c1，c2
            }
        }
        break;
    }
    return;
}

int main()
{
    // 初始化数据信息
    int data_num = 5, dimension = 2;
    vector<point> point_list, c1_data, c2_data;
    tmp_data(point_list);
    // random(point_list, data_num, dimension);
    vector<vector<int>> sorted_index = sort_data(point_list);

    // 划分原始数据，用于秘密共享
    divide_data(point_list, 10000, c1_data, c2_data);

    // 构造比较器
    Comparator *comparator = generate_comparator(false);

    // 初始化云
    cloud_one c1(c1_data, data_num, dimension, comparator, 5);
    cloud_two c2(c2_data, data_num, dimension, comparator, sorted_index, 5);

    // 构造乘法所需beaver三元组
    generate_beaver_set(data_num, 100, c1.beaver_list, c2.beaver_list);

    // 构造kd tree
    generate_kd_tree(c1, c2);

    // 聚类过程
    // filtering(c1, c2);

    return 0;
}