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

vector<Ctxt> Comparator::encrypt_vector(vector<int> x)
{
	// if(x.size()%2==1)
	// 	x.push_back(0);//FIXME 如果x长度为奇数，末尾补0，不会影响求最大值，但是会影响求最小值
	// scale参数标识是否对输入x进行缩放
	Ctxt ctxt_x(m_pk);
	vector<Ctxt> result(x.size(), ctxt_x);						  // 密文数组结果
	const EncryptedArray &ea = m_context.getEA();				  // 获取加密数组
	long nslots = ea.size();									  //! 提取槽slots的数量
	unsigned long p = m_context.getP();							  // p的值
	unsigned long ord_p = m_context.getOrdP();					  // p的次数
	unsigned long numbers_size = nslots / m_expansionLen;		  //! 一个密文中包含的数字
	unsigned long occupied_slots = numbers_size * m_expansionLen; //! 编码数字所需的槽数

	// encoding base, ((p+1)/2)^d
	// if 2-variable comparison polynomial is used, it must be p^d
	unsigned long enc_base = (p + 1) >> 1; // 俺懂了，p=7，那么(p+1)>>1 = 4， 所以二次编码的时候，就是以4为底
	if (m_type == BI || m_type == TAN)
		enc_base = p;
	unsigned long digit_base = power_long(enc_base, m_slotDeg); // m_slotDeg 是初始参数d，这里d=2

	// 检查field_size^expansion_len的值在64比特内
	int space_bit_size = static_cast<int>(ceil(m_expansionLen * log2(digit_base)));

	unsigned long input_range = ULONG_MAX;
	if (space_bit_size < 64)
		// input_range = power_long(field_size, expansion_len);
		input_range = power_long(digit_base, m_expansionLen); // 计算比较数字的最大范围

	cout << "最大输入：" << input_range << endl;
	cout << "一个密文可以包含的数字：" << numbers_size << endl;

	//! 存放加解密的结果
	vector<ZZX> expected_result(occupied_slots);
	vector<ZZX> decrypted(occupied_slots);

	// 对文本和模式构建明文多项式
	vector<ZZX> pol_x(nslots);

	//! 输入的内容
	unsigned long input_x;
	ZZX pol_slot;

	// FIXME 事情呢是这么个事儿
	// 没法加密负数，那么，只好，保证输入大于0吧！
	// 但是有点搞的是，加密整数，计算后如果是负数，再比较就没问题

	// 开始对输入进行拆分
	for (int i = 0; i < x.size(); ++i)
	{
		input_x = x[i];
		if (input_x < 0)
		{
			cout << "无法加密负数" << endl;
			exit(0);
		}
		// input_x = scale ? int(abs(x[i]) / 100) : x[i]; // 用scale参数来控制缩放

		if (input_x > input_range)
		{
			cout << input_x << " 数据超过加密范围..." << endl;
			exit(0);
		}

		if (m_verbose)
		{
			cout << "Input " << i << endl;
			cout << input_x << endl;
		} // 判断是否输出信息

		for (int j = 0; j < numbers_size; ++j)
		{

			vector<long> decomp_int_x;
			vector<long> decomp_char;

			digit_decomp(decomp_int_x, input_x, digit_base, m_expansionLen); // 分解数字

			if (m_verbose && !j)
			{
				cout << "Input decomposition into digits" << endl;
				for (int k = 0; k < m_expansionLen; k++)
					cout << decomp_int_x[k] << endl;
			} // 输出分解结果

			for (int k = 0; k < m_expansionLen; ++k)
			{ // 对槽进行编码
				// 分解数字
				int_to_slot(pol_slot, decomp_int_x[k], enc_base);
				pol_x[j * m_expansionLen + k] = pol_slot;
			}
		} // encode-->packing-->encrypt

		if (m_verbose)
		{
			cout << "Input" << endl;
			for (int i = 0; i < nslots; i++)
			{
				printZZX(cout, pol_x[i], ord_p);
				cout << endl;
			}
		} // 输出槽编码的结果

		//! 加密
		// Ctxt ctxt_x(m_pk);
		ea.encrypt(result[i], m_pk, pol_x);
		// result.push_back(ctxt_x);

	} // 每个数字加密为一个ciphertext，packing的数字都一样
	// Ctxt ct_res(m_pk);
	// compare(ct_res, result[0],result[1]);
	// result[0].multiplyBy(result[1]);
	// print_decrypted(result[0]);

	return result;
}

Ctxt Comparator::max_variance(vector<Ctxt> variance)
{
	// 构造所需的1密文 checked
	unsigned long p = m_context.getP();	   // p的值
	unsigned long enc_base = (p + 1) >> 1; // 俺懂了，p=7，那么(p+1)>>1 = 4， 所以二次编码的时候，就是以4为底
	if (m_type == BI || m_type == TAN)
		enc_base = p;
	unsigned long digit_base = power_long(enc_base, m_slotDeg); // m_slotDeg 是初始参数d，这里d=2

	int plain_one = 0;
	for (int i = 0; i < m_expansionLen; i++)
	{
		plain_one += int(pow(digit_base, i));
	}

	Ctxt ctxt_one(m_pk);
	ctxt_one = encrypt_vector(vector<int>{plain_one})[0];

	Ctxt max_value(m_pk);
	Ctxt max_index(m_pk);
	Ctxt less_than(m_pk);
	Ctxt mid_res(m_pk);

	vector<int> plain_index(variance.size());
	for (int i = 0; i < variance.size(); i++)
		plain_index[i] = i;
	vector<Ctxt> index = encrypt_vector(plain_index);
	vector<Ctxt> value = variance;

	while (value.size() != 1)
	{
		vector<Ctxt> new_value(value.size() / 2, ctxt_one);
		vector<Ctxt> new_index(value.size() / 2, ctxt_one);

		for (long int i = 0; i < value.size() / 2; i++)
		{
			compare(less_than, value[2 * i], value[2 * i + 1]);
			// update max value
			mid_res = ctxt_one;
			mid_res -= less_than;
			mid_res.multiplyBy(value[2 * i]);
			max_value = mid_res;

			mid_res = less_than;
			mid_res.multiplyBy(value[2 * i + 1]);
			max_value += mid_res;

			new_value[i] = max_value;

			// update max index
			mid_res = ctxt_one;
			mid_res -= less_than;
			mid_res.multiplyBy(index[2 * i]);
			max_index = mid_res;

			mid_res = less_than;
			mid_res.multiplyBy(index[2 * i + 1]);
			max_index += mid_res;

			new_index[i] = max_index;
		}

		if (value.size() % 2 == 1)
		{
			new_value.push_back(value.back());
			new_index.push_back(index.back());
		}

		value = new_value;
		index = new_index;
	}

	return index.back();
}