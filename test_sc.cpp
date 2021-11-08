#include <iostream>
#include <time.h>
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>
#include "tools.h"
#include "comparator.h"

using namespace std;
using namespace NTL;
using namespace helib;
using namespace he_cmp;
// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
// B 7 2 300 170 3 10 y
int main(){
	CircuitType type = BI;
	bool verbose = false;//是否输出
	unsigned long p = 7;
	unsigned long d = 2;	// 改变这个值就可以把【比较范围】变得好大好大，但是用时会增加
	unsigned long m = 300;
	unsigned long nb_primes = 500; // 噪声溢出修改这里
	unsigned long c = 3;

	auto context = ContextBuilder<BGV>()
	        .m(m)
	        .p(p)
	        .r(1)
	        .bits(nb_primes)
	        .c(c)
	        .scale(6)
	        .build();
	unsigned long expansion_len = 3;

	SecKey secret_key(context);
	secret_key.GenSecKey();

	// 生成密钥转换矩阵
	if (expansion_len > 1)
	{
	    if (context.getZMStar().numOfGens() == 1)
	    {
	      	std::set<long> automVals;
	      	long e = 1;
	      	long ord = context.getZMStar().OrderOf(0);
	      	bool native = context.getZMStar().SameOrd(0);
	      	if(!native)
	        	automVals.insert(context.getZMStar().genToPow(0, -ord));
	      	while (e < expansion_len){
	        	long atm = context.getZMStar().genToPow(0, ord-e);
	        	//cout << "Automorphism " << -e << " is " << atm << endl;
	        	automVals.insert(atm);
	        	e <<=1;
	      	}
	      	addTheseMatrices(secret_key, automVals);
	    }
	    else
	    {
	      addSome1DMatrices(secret_key);
	    }
	  }	
  	if (d > 1)
    	addFrbMatrices(secret_key); // might be useful only when d > 1

    // 创建比较器
	Comparator comparator(context, type, d, expansion_len, secret_key, verbose);

	int runs = 1;
	vector<int>tmp{1,3,2,5};

	vector<Ctxt>ctxt_vec = comparator.encrypt_vector(tmp);
	comparator.max_variance(ctxt_vec);

	return 0;

}