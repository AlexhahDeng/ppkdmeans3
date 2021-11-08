#include "my_tools.h"

//! 自定义工具函数
CircuitType type = BI;
bool verbose = false;//是否输出
unsigned long p = 7;
unsigned long d = 2;	// 改变这个值就可以把【比较范围】变得好大好大，但是用时会增加
unsigned long m = 300;
unsigned long nb_primes = 600; // 噪声溢出修改这里
unsigned long c = 3;

unsigned long expansion_len = 2;

auto context = ContextBuilder<BGV>()// 初始化上下文
	.m(m)// 分圆环的规模
	.p(p)// 明文素数模数
	.r(1)// 默认的什么东西
	.bits(nb_primes)// 模数链中密文素数比特的数量
	.c(c)// 密钥转换矩阵的列数
	.scale(6)// 不知道啥东西
	.build();

//! tools for comparison
Comparator *generate_comparator(bool output)
{
	/**
	 * func
	 * 		生成比较器
	 * para
	 * 		output:是否输出报错信息
	 * return
	 * 		指向比较器的指针
	 */

	CircuitType type = BI; 		// 比较电路的类型
	bool verbose = output;		// 是否输出错误信息
	unsigned long d = 2;		//! field extension degree 域扩展次数

	unsigned long expansion_len = 2;		  //!! maximal number of digits in a number
	SecKey *secret_key = new SecKey(context); // Create a secret key associated with the context
	secret_key->GenSecKey();				  // 生成加密所需私钥

	// 生成密钥转换矩阵
	if (expansion_len > 1)
	{
		if (context.getZMStar().numOfGens() == 1)
		{
			std::set<long> automVals;
			long e = 1;
			long ord = context.getZMStar().OrderOf(0);
			bool native = context.getZMStar().SameOrd(0);

			if (!native)
				automVals.insert(context.getZMStar().genToPow(0, -ord));

			while (e < expansion_len)
			{
				long atm = context.getZMStar().genToPow(0, ord - e);
				// cout << "Automorphism " << -e << " is " << atm << endl;
				automVals.insert(atm);
				e <<= 1;
			}
			addTheseMatrices(*secret_key, automVals);
		}
		else
			addSome1DMatrices(*secret_key);
	}

	if (d > 1)
		addFrbMatrices(*secret_key); // might be useful only when d > 1，不知道有啥用暂且先放在这里

	//! 创建比较器，然后发送给Cloud1，将公钥发送给User
	Comparator *comparator = new Comparator(context, type, d, expansion_len, *secret_key, verbose);

	return comparator;
}
