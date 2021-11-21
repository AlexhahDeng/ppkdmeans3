#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
using namespace std;

struct point{
	int index;
	vector<int>data;
};

void read_data(int data_num, int dimension){
	vector<point>result(data_num);
	ifstream fp("data.csv");
	string line;

	getline(fp, line);// 跳过列名
	for(int i=0;i<data_num;i++){// 循环读取每行数据
		getline(fp, line);
		vector<int>line_data(5);
		string number;

		istringstream readstr(line);//string数据流化
		// 将每一行按照 ，分割
		getline(readstr, number, ',');// 跳过行名

		for(int j=0;j<dimension;j++){
			getline(readstr, number, ',');
			line_data[j]=atoi(number.c_str());
		}

//		for(auto n:line_data){
//			cout<<n<<"  ";
//		}
//		cout<<endl;
		result[i]={i, line_data};
	}
}

void normalization(vector<point>point_list, int range){
	// 首先聚合point中的data
	vector<vector<int>>data_list(point_list.size());
	for(int i=0;i<point_list.size();i++){
		data_list[i]=point_list[i].data;
	}

	// 计算每个维度的最大值和最小值
	vector<vector<int>>max_min(data_list[0].size(), vector<int>(2));
	
}

int main(){
	read_data(8192, 5);
	return 0;
}
