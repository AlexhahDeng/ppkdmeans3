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

		for(auto n:line_data){
			cout<<n<<"  ";
		}
		cout<<endl;
		result[i]={i, line_data};
	}
}

int main(){
	read_data(8192, 5);
	return 0;
}
