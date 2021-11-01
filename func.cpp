#include "func.h"

void read_data()
{
	// 29 attributes, 65554 data
	ifstream inFile("../../data/Reaction Network (Undirected).data", ios::in);
	string lineStr;
	vector<vector<string>> strArray;
	// float array[65554][29] = { 0 };
	vector<vector<float>> array;

	if (inFile.fail())
		cout << "read failure!" << endl;

	while (getline(inFile, lineStr))
	{
		stringstream ss(lineStr);
		string str;
		vector<float> data;

		getline(ss, str, ','); // leave out the first one
		while (getline(ss, str, ','))
		{
			data.push_back(atof(str.c_str()));
		}
		array.push_back(data);
	}
	cout << array[0][0] << endl;
	cout << "number of instances: " << array.size() << endl;
	cout << "number of attributes: " << array[0].size() << endl;
}

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
			int num = rand() % point_list[i].data[j]; //! FIXME 这里没有模数据范围，而是模原始数据，保证秘密共享拆分结果都是整数，不知道会不会有问题
			c1[i].data[j] = num;
			c2[i].data[j] -= num;
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
		int a = rand() % data_range;
		int b = rand() % data_range;
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
	srand(time(NULL));
	// FIXME: 如果由cloud来随机选择初始簇中心，那不就暴露最初的中心了？虽然没有暴露值
	// 俺明白了，没关系！可以视为由用户初始化的啦！
	vector<int> ran_index(c1.k);
	int count = 0;

	while (count < c1.k)
	{
		int k_index = rand() % c1.data_num, isExist = 0;
		int i = 0;
		while (i < count)
		{
			if (ran_index[i] == k_index)
			{
				isExist = 1;
				break;
			} // 随机下标已经取过了，不能要
			i++;
		}
		if (isExist)
			continue;
		ran_index[count++] = k_index;
	}

	for (auto k_index : ran_index)
	{
		c1.clu_cen.push_back(c1.point_list[k_index].data);
		c2.clu_cen.push_back(c2.point_list[k_index].data);
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

	srand(time(NULL));
	for (int i = 0; i < c1.k; i++)
	{
		candidate1[i] = rand() % 2;
		candidate2[i] = 1 - candidate1[i];
	} //一开始候选集全部成立

	c1.kd_tree[0].candidate_k = candidate1;
	c2.kd_tree[0].candidate_k = candidate2;

	c1.clu_point_num = candidate1;
	c2.clu_point_num = candidate2;
}

vector<vector<int>> mul_in_vec(cloud_one &c1, cloud_two &c2, vector<int>& v1, vector<int>& v2)
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
	vector<int> r1 = c1.calculate_po_num(ef1);
	vector<int> r2 = c2.calculate_po_num(ef1);

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
	vector<vector<int>>vec{c1.clu_point_num, c2.clu_point_num};
	vector<int>r1(c1.k+1), r2(c2.k+1);

	for(int i = 0; i < c1.k; i++){
		vec[0][i]=0;
		vec[1][i]=1;
		vector<vector<int>>result(vec);

		while(true){
			vector<vector<int>>aha = mul_in_vec(c1, c2, result[0], result[1]);
			result = aha;
			if(aha[0].size() <= 1)
				break;
		}
		r1[i] = result[0][0];
		r2[i] = result[1][0];

		vec[0][i] = c1.clu_point_num[i];
		vec[1][i] = c2.clu_point_num[i];
	}

	// 全部相乘
	vector<vector<int>>result = vec;
	while(true){
		vector<vector<int>>aha = mul_in_vec(c1, c2, result[0], result[1]);
		result = aha;
		if(aha[0].size() <= 1)
			break;
	}
	r1[c1.k] = result[0][0];
	r2[c2.k] = result[1][0];
}