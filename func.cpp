#include "func.h"

void read_data() {
	// 29 attributes, 65554 data
	ifstream inFile("../../data/Reaction Network (Undirected).data", ios::in);
	string lineStr;
	vector<vector<string>> strArray;
	//float array[65554][29] = { 0 };
	vector<vector<float>>array;

	if (inFile.fail())
		cout << "read failure!" << endl;

	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		string str;
		vector<float>data;

		getline(ss, str, ','); // leave out the first one
		while (getline(ss, str, ',')) {
			data.push_back(atof(str.c_str()));
		}
		array.push_back(data);
	}
	cout << array[0][0] << endl;
	cout << "number of instances: " << array.size() << endl;
	cout << "number of attributes: " << array[0].size() << endl;
}

void quick_sort(vector<point>& s, int dimension, int l, int r) {
	/**
	* func: 对同一维度的数据进行排序
    * *checked
	*/
	if (l < r) {
		int i = l, j = r;
		point x = s[l];
		while (i < j) {
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

void random(vector<point>& point_list, int n, int m) {
	// number of data, number of dimension
	int l = 0, r = 1000;        // number range
	srand(time(NULL));

	for (int i = 0; i < n; i++) {

		vector<int>curr_data;
		for (int j = 0; j < m; j++)
			curr_data.push_back(rand() % (r - l + 1) + l);
		point_list.push_back({ i, curr_data });
	}
}

vector<vector<int>> sort_data(vector<point>point_list){
	/**
	 * func: sort data of every dimension
	 */
	int dimension = point_list[0].data.size();
	vector<vector<int>>sort_list;
	for(int i = 0; i < dimension; ++i){
		vector<int>res(point_list.size());

		quick_sort(point_list, i, 0, point_list.size() - 1);

		for(int j = 0; j < point_list.size(); ++j)
			res[j] = point_list[j].index;

		sort_list.push_back(res);
	}
	return sort_list;
}

void divide_data(vector<point>point_list, int data_range, vector<point>&c1, vector<point>&c2){
	/**
	 * func: 将数据随机拆分成两份
     * *checked
	 */
	c1 = point_list;
	c2 = point_list;
	srand(time(NULL));

	for(int i = 0; i < point_list.size(); ++i){
		for(int j = 0; j < point_list[0].data.size(); ++j){
			int num = rand() % point_list[i].data[j]; //! FIXME 这里没有模数据范围，而是模原始数据，保证秘密共享拆分结果都是整数，不知道会不会有问题
			c1[i].data[j] = num;
			c2[i].data[j] -= num;
		}
	}

}

void generate_beaver_set(int n, int data_range, vector<vector<int>>&c1_list, vector<vector<int>>&c2_list){
	/**
	 * func: 生成乘法所需beaver三元组
	 * note: c1_list, c2_list都是在云中生成好了的
	 * *checked
	 */
	srand(time(NULL));
	for(int i = 0; i < n; ++i){
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

void secret_share_N(vector<int>N, vector<int>&N1, vector<int>&N2){
	/**
	 * func: 将由01构成的N数组，划分为由01构成的N1，N2两个点集
     * *checked
	 */
	srand(time(NULL));
	for(int i = 0; i < N.size(); ++i){
		N1[i] = rand() % 2;
		N2[i] = N1[i] ^ N[i];
	}
}
