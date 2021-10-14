/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
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

// the main function that takes 7 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// argv[1] - circuit type (U, B or T)
// argv[2] - the plaintext modulus
// argv[3] - the dimension of a vector space over a finite field
// argv[4] - the order of the cyclotomic ring
// argv[5] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[6] - the length of vectors to be compared
// argv[7] - the number of experiment repetitions
// argv[8] - print debug info (y/n)

// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
int main(int argc, char *argv[]) {
  if(argc < 9)
  {
   throw invalid_argument("There should be exactly 8 arguments\n");
  }


  CircuitType type = UNI;
  // 判断使用的类型
  if (!strcmp(argv[1], "B"))
  {
    type = BI;
  }else if (!strcmp(argv[1], "T"))
  {
    type = TAN;
  }else if (!strcmp(argv[1], "U"))
  {
    type = UNI;
  }else
  {
    throw invalid_argument("Choose a valid circuit type (U for univariate, B for bivariate and T for Tan et al.\n");
  }

  // 是否输出全部信息，比如中间结果啥的
  bool verbose = false;
  if (!strcmp(argv[8], "y"))
    verbose = true;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus--明文模数
  unsigned long p = atol(argv[2]);
  // Field extension degree--域扩展次数
  unsigned long d = atol(argv[3]);
  // Cyclotomic polynomial - defines phi(m)--分圆多项式次数
  unsigned long m = atol(argv[4]);
  // Number of ciphertext prime bits in the modulus chain--模数链中密文素数比特的数量
  unsigned long nb_primes = atol(argv[5]);
  // Number of columns of Key-Switching matix (default = 2 or 3)--密钥转换矩阵的列数
  unsigned long c = 3;
  // 初始化上下文对象
  cout << "Initialising context object..." << endl;
  // Intialise context
  auto context = ContextBuilder<BGV>()
            .m(m)
            .p(p)
            .r(1)
            .bits(nb_primes)
            .c(c)
            .scale(6)
            .build();

  // Print the security level--输出安全等级
  // cout << "Q size: " << context.logOfProduct(context.getCtxtPrimes())/log(2.0) << endl;
  // cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
  // cout << "Security: " << context.securityLevel() << endl;

  // Print the context--输出上下文
  // context.getZMStar().printout();
  // cout << endl;

  //maximal number of digits in a number --一个数字中最大的digits数
  unsigned long expansion_len = atol(argv[6]);

  // Secret key management-- 生成私钥
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();

  // 生成密钥转换矩阵
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
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
    addFrbMatrices(secret_key); //might be useful only when d > 1

  // create Comparator (initialize after buildModChain)--创建比较器
  Comparator comparator(context, type, d, expansion_len, secret_key, verbose);

  //repeat experiments several times--重复几次实验
  int runs = atoi(argv[7]);
  
  //test comparison circuit -- 调用test_compare函数
  comparator.test_compare();

  //printAllTimers(cout);

  return 0;
}
