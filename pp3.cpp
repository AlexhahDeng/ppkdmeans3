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
#include <time.h>
#include "my_tools.h"

using namespace std;

/**
 * func： 构造kd tree
 */

void filtering2(cloud_one &c1, cloud_two &c2)
{

    ini_clu_cen(c1, c2);                           // 初始化簇中心（同时加密存放在c1中），存放在cloud中
    ini_candidate_k(c1, c2);                       // 初始化候选集和簇中点数，存放在kd_node中
    Ctxt ctxt_one = c1.comparator->gen_ctxt_one(); //目前没什么用，但是懒得改函数了-->用时将近0.001s
    Ctxt ctxt_zero = ctxt_one;                     // 构造密文0
    ctxt_zero *= 0l;
    Ctxt less_than = ctxt_one;

    int iter = 0;
    while (iter < 5)
    {
        cout << "**********iteration" << iter << "**********" << endl;
        clock_t sta, end;
        sta = clock();

        mul_clu_point_num(c1, c2);                        //连乘α
        c1.kd_tree[0].candidate_k = vector<int>(c1.k, 1); //心中随机即可
        c2.kd_tree[0].candidate_k = vector<int>(c1.k, 0); // 初始化root候选中心

        vector<vector<long int>> knum_square = cal_knum_square(c1, c2); // 计算簇包含点数的平方-->checked

        vector<int> new_clu_point_num1(c1.k, 0), new_clu_point_num2(new_clu_point_num1);
        vector<vector<int>> new_clu_cen1(c1.k, vector<int>(c1.dimension, 0)), new_clu_cen2(new_clu_cen1); //? 这里可以把密文结果累加，存储到一起，迭代完再进行ss划分成两份

        // 开始迭代

        for (int node_index = 0; node_index < c1.kd_tree.size(); node_index++)
        {
            // cout << "*******************************" << endl
            //      << "正在聚类node " << node_index << endl;

            if (c1.kd_tree[node_index].isClustered)
            {
                if (2 * node_index + 1 < c1.kd_tree.size())
                    c1.kd_tree[2 * node_index + 1].isClustered = true;
                if (2 * node_index + 2 < c2.kd_tree.size())
                    c1.kd_tree[2 * node_index + 2].isClustered = true;
                continue;
            } // 已被划分则不再参与划分

            int clo_k_index = 0;
            vector<long int> dist = cal_dist(c1, c2, node_index, 0);                   //计算node到每个簇的距离，非candidate的簇距离为0
            vector<Ctxt> ctxt_dist = c1.comparator->encrypt_dist(dist, clo_k_index);   //加密距离
            vector<Ctxt> min_dist_mark = c1.comparator->min_dist(ctxt_dist, ctxt_one); //求最近簇的标识

            if (c1.kd_tree[node_index].node_point_num > 3)
            { // not leaf node
                vector<Ctxt> ctxt_clo_k(c1.dimension, ctxt_zero);
                for (int j = 0; j < c1.k; j++)
                {
                    for (int k = 0; k < c1.dimension; k++)
                    {
                        Ctxt ctxt_value = min_dist_mark[j]; // 取当前k对应的min dist计算结果
                        ctxt_value *= c1.ctxt_clu_cen[j][k];
                        ctxt_clo_k[k] += ctxt_value;
                    }
                } // 计算z*-->距离最近的簇

                // 最近中心转ss
                // vector<int> clo_k1(c1.dimension), clo_k2(clo_k1);
                // c1.comparator->he_to_ss(ctxt_clo_k, clo_k1, clo_k2);

                // 开始遍历所有簇
                for (int k_index = 0; k_index < c1.k; k_index++)
                {
                    //cout << "正在判断簇 " << k_index << endl;
                    vector<long int> v1(c1.dimension), v2(v1);

                    for (int j = 0; j < c1.dimension; j++)
                    {
                        c1.comparator->compare(less_than, c1.ctxt_clu_cen[k_index][j], ctxt_clo_k[j]); // zij < kj ?
                        long com_res = c1.comparator->dec_compare_res(less_than);

                        int mid_res = com_res * c1.kd_tree[node_index].ptxt_node_min[j] + (1 - com_res) * c1.kd_tree[node_index].ptxt_node_max[j];
                        v1[j] = rand() % (mid_res + 1);
                        v2[j] = mid_res - v1[j];
                    }

                    // 计算v*|z|
                    vector<vector<long int>> v_mul_znum = cal_v_znum(c1, c2, v1, v2, k_index); //->checked
                    vector<vector<long int>> v_mul_clok = cal_v_znum(c1, c2, v1, v2, clo_k_index);

                    // 计算 Σ（z-v|z|)**2-->res=2*dimension
                    vector<long int> k_min_znum = cal_vz_sqaure(c1, c2, v_mul_znum, k_index); //->checked
                    vector<long int> clo_min_znum = cal_vz_sqaure(c1, c2, v_mul_clok, clo_k_index);

                    // 计算距离-->checked
                    vector<long int> distkv = mul_two(c1, c2, k_min_znum, vector<long int>{knum_square[0][k_index], knum_square[1][k_index]});
                    vector<long int> distclo_kv = mul_two(c1, c2, clo_min_znum, vector<long int>{knum_square[0][clo_k_index], knum_square[1][clo_k_index]});

                    // 加密距离并进行比较
                    vector<Ctxt> enc_dist = c1.comparator->encrypt_dist_vk(vector<long int>{distkv[0] + distkv[1], distclo_kv[0] + distclo_kv[1]});
                    c1.comparator->compare(less_than, enc_dist[1], enc_dist[0]); // dist[z*,v]<dist[z,v]?
                    long com_res = c1.comparator->dec_compare_res(less_than);

                    c1.kd_tree[node_index].candidate_k[k_index] *= (1 - com_res);
                    if (node_index * 2 + 1 < c1.kd_tree.size())
                        c1.kd_tree[node_index * 2 + 1].candidate_k[k_index] = c1.kd_tree[node_index].candidate_k[k_index];
                    if (node_index * 2 + 2 < c1.kd_tree.size())
                        c1.kd_tree[node_index * 2 + 2].candidate_k[k_index] = c1.kd_tree[node_index].candidate_k[k_index];
                }

                int tot_can_k = 0;
                for (int i = 0; i < c1.k; i++)
                    tot_can_k += c1.kd_tree[node_index].candidate_k[i];

                if (tot_can_k <= 1)
                {
                    cout << "node " << node_index << " is clustered" << endl;

                    vector<int> ptxt_min_mark1(c1.k, 0), ptxt_min_mark2(c1.k, 0); // 秘密共享最近簇标识
                    c1.comparator->dist_mark_to_ss(min_dist_mark, ptxt_min_mark1, ptxt_min_mark2);

                    vector<vector<int>> e1, f1, e2, f2;
                    c1.update_clu_cen_ef(node_index, ptxt_min_mark1, e1, f1);
                    c2.update_clu_cen_ef(node_index, ptxt_min_mark2, e2, f2);

                    for (int i = 0; i < e1.size(); i++)
                    {
                        for (int j = 0; j < e1[0].size(); j++)
                        {
                            e1[i][j] += e2[i][j];
                            f1[i][j] += f2[i][j];
                        }
                    }

                    vector<vector<int>> r1 = c1.update_clu_cen_final(e1, f1);
                    vector<vector<int>> r2 = c2.update_clu_cen_final(e1, f1); //size=(dimension+1)*k

                    for (int i = 0; i < c1.k; i++)
                    {
                        for (int j = 0; j < r1.size() - 1; j++)
                        {
                            new_clu_cen1[i][j] = r1[j][i];
                            new_clu_cen2[i][j] = r2[j][i];
                        }
                        new_clu_point_num1[i] += c1.kd_tree[node_index].node_point_num * (ptxt_min_mark1[i]+ptxt_min_mark2[i]);
                        new_clu_point_num2[i] += 0;
                    } // i 对应k，j对应维度数据

                    // 标识节点信息
                    c1.kd_tree[node_index].isClustered = true;
                    if ((node_index * 2 + 1) < c1.kd_tree.size())
                        c1.kd_tree[node_index * 2 + 1].isClustered = true;
                    if ((node_index * 2 + 2) < c1.kd_tree.size())
                        c1.kd_tree[node_index * 2 + 2].isClustered = true;

                } // FIXME好像也可以切成秘密共享，麻了，再说吧
            }
            else
            { // leaf node
                cout << "leaf node " << node_index << " is clustered" << endl;

                vector<int> ptxt_min_mark1(c1.k, 0), ptxt_min_mark2(c1.k, 0); // 秘密共享最近簇标识
                c1.comparator->dist_mark_to_ss(min_dist_mark, ptxt_min_mark1, ptxt_min_mark2);

                vector<vector<int>> e1, f1, e2, f2;
                c1.update_clu_cen_ef(node_index, ptxt_min_mark1, e1, f1);
                c2.update_clu_cen_ef(node_index, ptxt_min_mark2, e2, f2);

                for (int i = 0; i < e1.size(); i++)
                {
                    for (int j = 0; j < e1[0].size(); j++)
                    {
                        e1[i][j] += e2[i][j];
                        f1[i][j] += f2[i][j];
                    }
                }

                vector<vector<int>> r1 = c1.update_clu_cen_final(e1, f1);
                vector<vector<int>> r2 = c2.update_clu_cen_final(e1, f1); //size=(dimension+1)*k

                for (int i = 0; i < c1.k; i++)
                {
                    for (int j = 0; j < r1.size() - 1; j++)
                    {
                        new_clu_cen1[i][j] += r1[j][i];
                        new_clu_cen2[i][j] += r2[j][i];
                    }
                    new_clu_point_num1[i] += c1.kd_tree[node_index].node_point_num * (ptxt_min_mark1[i]+ptxt_min_mark2[i]);
                    new_clu_point_num2[i] += 0;
                } // i 对应k，j对应维度数据

                c1.kd_tree[node_index].isClustered = true;
                if ((node_index * 2 + 1) < c1.kd_tree.size())
                    c1.kd_tree[node_index * 2 + 1].isClustered = true;
                if ((node_index * 2 + 2) < c1.kd_tree.size())
                    c1.kd_tree[node_index * 2 + 2].isClustered = true;
            }
        }
        end = clock();
        cout << "tot time " << (double)(end - sta) / CLOCKS_PER_SEC << "s" << endl;

        c1.clu_cen = new_clu_cen1;
        c2.clu_cen = new_clu_cen2;
        c1.clu_point_num = new_clu_point_num1;
        c2.clu_point_num = new_clu_point_num2;

        for(int i=0;i<c1.kd_tree.size();i++){
            c1.kd_tree[i].isClustered = false;
        }

        iter++;
    }
    return;
}
void generate_kd_tree2(cloud_one &c1, cloud_two &c2)
{

    c2.add_new_node(vector<int>(c2.data_num, 1), c1.data_num, c1.kd_tree, c2.kd_tree);
    int node_index = 0;

    while (node_index < c1.kd_tree.size())
    {
        cout << "生成kd tree 节点 " << node_index << endl;

        int point_num = c1.kd_tree[node_index].node_point_num; // 获取当前tree node包含的点个数

        //* 计算xi*Ni
        // 1. 初始化e1，e2 和 f1，f2，二维数组，data_num*dimension
        vector<vector<int>> e1(c1.data_num, vector<int>(c1.dimension, 0)), e2(e1);
        vector<vector<int>> f1(e1), f2(e1);

        // 2. 计算xi*Ni的中间结果 e1,f1 and e2, f2
        c1.calculate_ef_xN(node_index, e1, f1);
        c2.calculate_ef_xN(node_index, e2, f2);

        // 3. 计算合并结果 e=e1+e2, f=f1+f2,那么就用e1和f1暂存中间结果e，f吧！
        for (int i = 0; i < c1.data_num; ++i)
        {
            for (int j = 0; j < c1.dimension; ++j)
            {
                e1[i][j] += e2[i][j];
                f1[i][j] += f2[i][j];
            }
        }

        // 4. xN1 = f1*a1+e1*b1+c1  xN2=e2*f2+f2*a2+e2*b2+c2 --> checked
        vector<vector<int>> xN1 = c1.calculate_xi_Ni(e1, f1);
        vector<vector<int>> xN2 = c2.calculate_xi_Ni(e1, f1);

        //* 计算平均值(ΣxN1)/n (ΣxN2)/n
        vector<int> avg_xN1(c1.dimension, 0), avg_xN2(c1.dimension, 0);
        vector<int> ptxt_xN(avg_xN1);

        //-->求和，取平均
        for (int i = 0; i < c1.dimension; i++)
        {
            for (int j = 0; j < c1.data_num; j++)
            {
                avg_xN1[i] += xN1[j][i];
                avg_xN2[i] += xN2[j][i];
            } // 也可以边累加，边做除法，误差比较大，但是中间结果不会溢出

            //! 将(Σxi*Ni)补充到kd_node中
            c1.kd_tree[node_index].node_sum_x[i] = avg_xN1[i];
            c2.kd_tree[node_index].node_sum_x[i] = avg_xN2[i];

            ptxt_xN[i] = avg_xN1[i] + avg_xN2[i]; // 记录明文Σxi*Ni，不可以是均值！

            avg_xN1[i] /= point_num;
            avg_xN2[i] /= point_num;

        } // 果然最后做除法，误差比“边加边除”方法小到了小数级别

        c1.kd_tree[node_index].ctxt_node_sum = c1.comparator->encrypt_vector(ptxt_xN);
        //! 在节点中，补充节点密文和，方便后续计算-->有一定可能超过范围哦-->也不一定，就是说，看我把数据放大的范围吧

        if (point_num <= 3)
        {
            node_index++;
            continue;
        } // 叶子节点只记录ΣxiNi，不做其他计算

        // 计算avg*Ni--> avg(data_num*dimension)
        c1.calculate_avg_N(e1, f1, avg_xN1, node_index);
        c2.calculate_avg_N(e2, f2, avg_xN2, node_index);

        //-->合并中间结果
        for (int i = 0; i < c1.data_num; i++)
        {
            for (int j = 0; j < c1.dimension; j++)
            {
                e1[i][j] += e2[i][j];
                f1[i][j] += f2[i][j];
            }
        }

        //-->计算avg*Ni最终结果
        vector<vector<int>> avg_mul_N1 = c1.calculate_xi_Ni(e1, f1);
        vector<vector<int>> avg_mul_N2 = c2.calculate_xi_Ni(e1, f1);

        // 计算xNi - avg_mul_Ni --> 暂存到前者中
        for (int i = 0; i < c1.data_num; i++)
        {
            for (int j = 0; j < c1.dimension; j++)
            {
                xN1[i][j] -= avg_mul_N1[i][j];
                xN2[i][j] -= avg_mul_N2[i][j];
            }
        }

        // 计算 (xNi - avg_xNi)^2
        //--> 计算中间结果ef(借用现有函数)
        c1.calculate_ef_xN_square(xN1, e1, f1);
        c2.calculate_ef_xN_square(xN2, e2, f2);

        //--> 合并中间结果ef，暂存到e1,f1中
        for (int i = 0; i < c1.data_num; i++)
        {
            for (int j = 0; j < c1.dimension; j++)
            {
                e1[i][j] += e2[i][j];
                f1[i][j] += f2[i][j];
            }
        }

        // 计算方差最终--> checked
        vector<vector<int>> xN_sub_avg1 = c1.calculate_xi_Ni_square(e1, f1);
        vector<vector<int>> xN_sub_avg2 = c2.calculate_xi_Ni_square(e1, f1);

        //* 比较方差大小
        vector<int> var(c1.dimension, 0);
        for (int i = 0; i < c1.dimension; i++)
        {
            float curr = 0;
            for (int j = 0; j < c1.data_num; j++)
                curr += float(xN_sub_avg1[j][i] + xN_sub_avg2[j][i]) / point_num; // 直接合并秘密共享的中间值了（原本应当分开加密再合并的
            var[i] = curr;
        } // 看了下，如果合并再除法，中间结果快上亿了，估计会有运算错误，因此一边除法一边合并

        // 加密，比较大小
        vector<Ctxt> enc_var = c1.comparator->encrypt_vector(var); // 加密方差
        Ctxt ctxt_one = c1.comparator->gen_ctxt_one();             // 生成密文全1

        Ctxt max_var_index = c1.comparator->max_variance(enc_var, ctxt_one); // 比较获取最大方差的维度
        // cout<<c1.comparator->decrypt_index(max_var_index);                            // 解密最大方差维度的下标

        vector<int> N(c1.data_num);
        for (int i = 0; i < c1.data_num; i++)
        {
            N[i] = c1.kd_tree[node_index].N[i] + c2.kd_tree[node_index].N[i];
        }

        c2.divide_data_set(max_var_index, c1, N, point_num, node_index); // 根据最大方差的维度划分数据

        node_index++;
    }

    return;
}
int main()
{
    setbuf(stdout, 0);
    srand(time(NULL));
    // 初始化数据信息
    int data_num = 100, dimension = 5;
    vector<point> point_list, c1_data, c2_data;

    // 读取数据
    point_list = read_data(data_num, dimension);

    // user：对数据做排序
    vector<vector<int>> sorted_index = sort_data(point_list);

    // user：划分原始数据，用于秘密共享
    divide_data(point_list, 10000, c1_data, c2_data);

    // 构造比较器
    Comparator *comparator = generate_comparator(false);

    // 初始化云
    cloud_one c1(c1_data, data_num, dimension, comparator, 3);
    cloud_two c2(c2_data, data_num, dimension, comparator, sorted_index, 3);

    // 构造乘法所需beaver三元组
    generate_beaver_set(data_num, 10, c1.beaver_list, c2.beaver_list);

    // 构造kd tree
    generate_kd_tree2(c1, c2);

    // 聚类过程
    filtering2(c1, c2);

    return 0;
}
