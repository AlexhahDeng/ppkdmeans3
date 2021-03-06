#include "func.h"

vector<point> read_data(int data_num, int dimension)
{
	vector<point> result(data_num);
	ifstream fp("data.csv");
	string line;

	// getline(fp, line); // 跳过列名-->新数据没有列名
	for (int i = 0; i < data_num; i++)
	{ // 循环读取每行数据

		getline(fp, line);
		vector<int> line_data(5);
		string number;

		istringstream readstr(line); // string数据流化
		// 将每一行按照 "," 分割
		// getline(readstr, number, ','); // 跳过行名-->新数据没有行名

		for (int j = 0; j < dimension; j++)
		{
			getline(readstr, number, ',');
			line_data[j] = atoi(number.c_str());
		}

		// for(auto n:line_data){
		// 	cout<<n<<"  ";
		// }
		// cout<<endl;
		result[i] = {i, line_data};
	}
	fp.close();
	return result;
}

// void write_file(vector<point>p_list, vector<point>c1, vector<point>c2){

// }

void quick_sort(vector<point> &s, int dimension, int l, int r)
{
	/**
	 * func: 对同一维度的数据进行排序
	 * *checked
	 */
	if (l < r)
	{
		int i = l, j = r;
		point x = s[l];
		while (i < j)
		{
			while (i < j && s[j].data[dimension] >= x.data[dimension])
				j--;

			if (i < j)
				s[i++] = s[j];

			while (i < j && s[i].data[dimension] < x.data[dimension])
				i++;

			if (i < j)
				s[j--] = s[i];
		}
		s[i] = x;
		quick_sort(s, dimension, l, i - 1);
		quick_sort(s, dimension, i + 1, r);
	}
}

void random(vector<point> &point_list, int n, int m)
{
	// number of data, number of dimension
	int l = 0, r = 1000; // number range
	srand(time(NULL));

	for (int i = 0; i < n; i++)
	{

		vector<int> curr_data;
		for (int j = 0; j < m; j++)
			curr_data.push_back(rand() % (r - l + 1) + l);
		point_list.push_back({i, curr_data});
	}
}

vector<vector<int>> sort_data(vector<point> point_list)
{
	/**
	 * func: sort data of every dimension
	 */
	int dimension = point_list[0].data.size();
	vector<vector<int>> sort_list;
	for (int i = 0; i < dimension; ++i)
	{
		vector<int> res(point_list.size());

		quick_sort(point_list, i, 0, point_list.size() - 1);

		for (int j = 0; j < point_list.size(); ++j)
			res[j] = point_list[j].index;

		sort_list.push_back(res);
	}
	return sort_list;
}

void divide_data(vector<point> point_list, int data_range, vector<point> &c1, vector<point> &c2)
{
	/**
	 * func: 将数据随机拆分成两份
	 * *checked
	 */
	c1 = point_list;
	c2 = point_list;
	srand(time(NULL));

	for (int i = 0; i < point_list.size(); ++i)
	{
		for (int j = 0; j < point_list[0].data.size(); ++j)
		{
			if (point_list[i].data[j] == 0) // 如果数据为0， 那就都为0
				continue;
			else
			{
				int num = rand() % point_list[i].data[j]; // 模的是原始数据哦，保证秘密共享拆分都是正整数
				c1[i].data[j] = num;
				c2[i].data[j] -= num;
			}
		}
	}
}

void generate_beaver_set(int n, int data_range, vector<vector<int>> &c1_list, vector<vector<int>> &c2_list)
{
	/**
	 * func: 生成乘法所需beaver三元组
	 * note: c1_list, c2_list都是在云中生成好了的
	 * *checked
	 */
	srand(time(NULL));
	for (int i = 0; i < n; ++i)
	{
		int a = rand() % data_range + 1;
		int b = rand() % data_range + 1;
		int c = a * b;

		int c1_a = rand() % a;
		int c2_a = a - c1_a;

		int c1_b = rand() % b;
		int c2_b = b - c1_b;

		int c1_c = rand() % c;
		int c2_c = c - c1_c;

		c1_list[i] = vector<int>{c1_a, c1_b, c1_c};
		c2_list[i] = vector<int>{c2_a, c2_b, c2_c};
	}
}

void tmp_data(vector<point> &point_list)
{
	/**
	 * func: 固定化测试数据
	 */
	point_list.push_back({0, vector<int>{123, 43}});
	point_list.push_back({1, vector<int>{76, 987}});
	point_list.push_back({2, vector<int>{234, 76}});
	point_list.push_back({3, vector<int>{67, 567}});
	point_list.push_back({4, vector<int>{59, 473}});
}

// func: 用户初始化簇中心，确保点不重复--checked
void ini_clu_cen(cloud_one &c1, cloud_two &c2)
{
	// 如果由cloud来随机选择初始簇中心，那不就暴露最初的中心了？虽然没有暴露值
	// SOLUTION: 俺明白了，没关系！可以视为由用户初始化的啦！
	//	vector<int> ran_index(c1.k);
	//	int count = 0;
	//
	//	while (count < c1.k)
	//	{
	//		int k_index = rand() % c1.data_num, isExist = 0;
	//		int i = 0;
	//		while (i < count)
	//		{
	//			if (ran_index[i] == k_index)
	//			{
	//				isExist = 1;
	//				break;
	//			} // 随机下标已经取过了，不能要
	//			i++;
	//		} // 检测随机选择的簇中心是否已经被选择过
	//		if (isExist)
	//			continue;
	//		ran_index[count++] = k_index;
	//	}
//	vector<int> ran_index{1, 2, 3};

    vector<int> ran_index(c1.k);
    for(int i=0;i<ran_index.size();i++){
        ran_index[i] = rand() % c1.data_num;
    }
	for (auto k_index : ran_index)
	{
		c1.clu_cen.push_back(c1.point_list[k_index].data);
		c2.clu_cen.push_back(c2.point_list[k_index].data);

		vector<int> k_cen(c1.dimension);
		for (int i = 0; i < c1.dimension; i++)
			k_cen[i] = c1.point_list[k_index].data[i] + c2.point_list[k_index].data[i];

		c1.ctxt_clu_cen.push_back(c1.comparator->encrypt_vector(k_cen));
		// 这里把簇中心加密，存入c1，方便聚类过程中的计算
	}

	return;
}

// func: 用户初始化候选中心和簇大小给root--checked
void ini_candidate_k(cloud_one &c1, cloud_two &c2)
{
	// 初始每个簇中包含的点数为1，这样不影响公式的正常计算！
	// 刚好vector的长度和candidate set一样，那就一家人不说两家话，弄成一样的呗
	vector<int> candidate1(c1.k);
	vector<int> candidate2(c2.k);

	for (int i = 0; i < c1.k; i++)
	{
		candidate1[i] = 1;
		candidate2[i] = 0;
	} //一开始候选集全部成立-->心中随机

	c1.kd_tree[0].candidate_k = candidate1;
	c2.kd_tree[0].candidate_k = candidate2;

	c1.kd_tree[1].candidate_k = candidate1;
	c2.kd_tree[1].candidate_k = candidate2;

	c1.kd_tree[2].candidate_k = candidate1;
	c2.kd_tree[2].candidate_k = candidate2;

	c1.clu_point_num = candidate1;
	c2.clu_point_num = candidate2;
}

vector<vector<int>> mul_in_vec(cloud_one &c1, cloud_two &c2, vector<int> &v1, vector<int> &v2)
{
	// 向量v中所有奇偶元素相乘，也就是01，23，45...
	// 最后剩下的元素直接添加到result末尾
	// rol0: result1   rol1: result2  size = 2*len
	// *checked
	vector<vector<int>> ef1 = c1.calculate_po_num_ef(v1);
	vector<vector<int>> ef2 = c2.calculate_po_num_ef(v2);

	for (int i = 0; i < ef1.size(); i++)
	{
		ef1[i][0] += ef2[i][0];
		ef1[i][1] += ef2[i][1];
	} // 合并e，f，暂存在1中

	// 计算结果
	vector<int> r1 = c1.calculate_mul_final(ef1);
	vector<int> r2 = c2.calculate_mul_final(ef1);

	if (v1.size() % 2 != 0)
	{
		r1.push_back(v1.back());
		r2.push_back(v2.back());
	} // 如果维度是奇数，把最后一个没乘的补上

	return vector<vector<int>>{r1, r2};
}

void mul_clu_point_num(cloud_one &c1, cloud_two &c2)
{
	// TODO 烦死啦！只能两个两个乘，二分把！
	vector<vector<int>> vec{c1.clu_point_num, c2.clu_point_num};
	vector<int> r1(c1.k + 1), r2(c2.k + 1);

	for (int i = 0; i < c1.k; i++)
	{
		vec[0][i] = 0;
		vec[1][i] = 1;
		vector<vector<int>> result(vec);

		while (true)
		{
			vector<vector<int>> aha = mul_in_vec(c1, c2, result[0], result[1]);
			result = aha;
			if (aha[0].size() <= 1)
				break;
		}
		r1[i] = result[0][0];
		r2[i] = result[1][0];

		vec[0][i] = c1.clu_point_num[i];
		vec[1][i] = c2.clu_point_num[i];
	}

	// 全部相乘
	vector<vector<int>> result = vec;
	while (true)
	{
		vector<vector<int>> aha = mul_in_vec(c1, c2, result[0], result[1]);
		result = aha;
		if (aha[0].size() <= 1)
			break;
	}
	r1[c1.k] = result[0][0];
	r2[c2.k] = result[1][0];

	c1.mul_point_num = r1;
	c2.mul_point_num = r2;
	return;
}

// 直接返回合并秘密共享的结果，以及和n乘的结果
vector<long int> cal_dist(cloud_one &c1, cloud_two &c2, int node_index, int tot_can_num)
{

	vector<long int> result(c1.k, 0);
	//	clock_t start, end;
	//	start = clock();
	for (int i = 0; i < c1.k; ++i)
	{
		vector<vector<long int>> ef1 = c1.calculate_dist_para_ef(node_index, i);
		vector<vector<long int>> ef2 = c2.calculate_dist_para_ef(node_index, i);

		for (int j = 0; j < ef1.size(); j++)
		{
			ef1[j][0] += ef2[j][0];
			ef1[j][1] += ef2[j][1];
		} // 合并e，f，暂存到ef1中

		// 获取中间参数的乘法结果  Σc^2, Σk^2, ci*ki, α^2, αj^2, ααj
		vector<long int> para1 = c1.calculate_dist_para(ef1);
		vector<long int> para2 = c2.calculate_dist_para(ef1); // 我好像找到bug了！！！

		// 计算distance中间结果ef
		ef1 = c1.calculate_dist_res_ef(para1);
		ef2 = c2.calculate_dist_res_ef(para2);

		// 合并e和f
		for (int j = 0; j < ef1.size(); j++)
		{
			ef1[j][0] += ef2[j][0];
			ef1[j][1] += ef2[j][1];
		}

		// 计算到簇k距离最终结果
		int n = c1.kd_tree[node_index].candidate_k[i] + c2.kd_tree[node_index].candidate_k[i]; // 是否已被prune
		tot_can_num += n;
		result[i] = (c1.calculate_dist_res(ef1) + c2.calculate_dist_res(ef1)) * n;

	} // 计算中心到第i个簇的距离

	return result;
}

vector<vector<long int>> cal_v_znum(cloud_one &c1, cloud_two &c2, vector<long int> v1, vector<long int> v2, int k_index)
{
	// 计算ef
	vector<vector<long int>> ef1 = c1.cal_vznum_ef(v1, k_index);
	vector<vector<long int>> ef2 = c2.cal_vznum_ef(v2, k_index);

	//合并ef到ef1中
	for (int i = 0; i < ef1.size(); i++)
	{
		ef1[i][0] += ef2[i][0];
		ef1[i][1] += ef2[i][1];
	}

	// 获取计算结果
	vector<long int> r1 = c1.cal_vznum_final(ef1);
	vector<long int> r2 = c2.cal_vznum_final(ef1);

	return vector<vector<long int>>{r1, r2};
}

vector<long int> cal_vz_sqaure(cloud_one &c1, cloud_two &c2, vector<vector<long int>> vz, int k_index)
{
	// vz--> 2xdimension
	vector<int> tmp1 = c1.clu_cen[k_index];
	vector<int> tmp2 = c2.clu_cen[k_index];

    vector<long int> r1(tmp1.size()), r2(r1);

	for (int i = 0; i < c1.dimension; i++)
	{
		r1[i] = tmp1[i] - vz[0][i];
		r2[i] = tmp2[i] - vz[1][i];
	} // z[i] - v*|z|[i]

	// 计算r[i]**2
	vector<vector<long int>> ef1 = c1.cal_vec_square(r1);
	vector<vector<long int>> ef2 = c2.cal_vec_square(r2);

	for (int i = 0; i < ef1.size(); i++)
	{
		ef1[i][0] += ef2[i][0];
		ef1[i][1] += ef2[i][1];
	}

	// 计算最终结果
	r1 = c1.cal_square_final(ef1);
	r2 = c2.cal_square_final(ef1);

	return vector<long int>{accumulate(r1.begin(), r1.end(), 0), accumulate(r2.begin(), r2.end(), 0)};
}
// 计算每个簇包含点数的平方，vec.size()=k
vector<vector<long int>> cal_knum_square(cloud_one &c1, cloud_two &c2)
{
	vector<int> knum1 = c1.clu_point_num;
	vector<int> knum2 = c2.clu_point_num;

	vector<vector<long int>> ef1 = c1.cal_vec_square(knum1);
	vector<vector<long int>> ef2 = c2.cal_vec_square(knum2);

	for (int i = 0; i < ef1.size(); i++)
	{
		ef1[i][0] += ef2[i][0];
		ef1[i][1] += ef2[i][1];
	}

	// 计算最终结果
	vector<long int> r1 = c1.cal_square_final(ef1);
	vector<long int> r2 = c2.cal_square_final(ef1);

	return vector<vector<long int>>{r1, r2};
}

vector<long int> mul_two(cloud_one &c1, cloud_two &c2, vector<long int> numa, vector<long int> numb)
{
	//计算两个数的乘积，numa和numb分别都是秘密共享值（size=2-->r0=ss1,r0=ss2

	// 计算ef
	long int e1 = numa[0] - c1.beaver_list[0][0]; // numa1-a
	long int f1 = numb[0] - c1.beaver_list[0][1]; // numb1-b

	long int e2 = numa[1] - c2.beaver_list[0][0]; // numa2-a
	long int f2 = numb[1] - c2.beaver_list[0][1]; // numa2-b

	long int e = e1 + e2;
	long int f = f1 + f2;

	long int r1 = e * c1.beaver_list[0][1] + f * c1.beaver_list[0][0] + c1.beaver_list[0][2];
	long int r2 = e * f + e * c2.beaver_list[0][1] + f * c2.beaver_list[0][0] + c2.beaver_list[0][2];

	return vector<long int>{r1, r2};
}
