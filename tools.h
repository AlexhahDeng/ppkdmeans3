#pragma once
/*
Auxiliary functions for integer encoding
存放一些工具函数
*/

#ifndef TOOLS_H
#define TOOLS_H
#include <iostream>
#include <cmath>
#include <vector>
#include <helib/helib.h>
#include <helib/Ctxt.h>
#include <helib/polyEval.h>
#include <helib/debugging.h>
#include <random>
#include <map> 
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <helib/Ptxt.h>
#include "comparator.h"

using namespace std;
using namespace helib;
using namespace NTL;
using namespace he_cmp;

//! defined by myself
Comparator* generate_comparator();

//! tools for comparison
inline int intlog(unsigned long base, unsigned long input)
{
	int res = max(static_cast<int>(floor(log2(input) / log2(base))), 0);
	if (power_long(base, res + 1) == input)
		return res + 1;
	return res;
}

void digit_decomp(vector<long>& decomp, unsigned long input, unsigned long base, int nslots);

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
void simplePolyEval(Ctxt& ret, const NTL::ZZX& poly, DynamicCtxtPowers& babyStep);

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
void PatersonStockmeyer(Ctxt& ret, const NTL::ZZX& poly, long k, long t, long delta, DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
void degPowerOfTwo(Ctxt& ret, const NTL::ZZX& poly, long k,
	DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

void recursivePolyEval(Ctxt& ret, const NTL::ZZX& poly, long k,
	DynamicCtxtPowers& babyStep, DynamicCtxtPowers& giantStep);

#endif // #ifndef TOOLS_H