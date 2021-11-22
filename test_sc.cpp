#include <iostream>
#include <time.h>
#include "my_tools.h"
#include <iomanip>
// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
// B 7 2 300 170 3 10 y
int main(){
	Comparator *comparator = generate_comparator(false);

	int runs = 1;

	vector<Ctxt>test();
	vector<long int>tmp{0,117648,49,255,5};
	Ctxt ctxt_one = comparator->gen_ctxt_one();

	vector<Ctxt>ctxt_vec = comparator->encrypt_vector(tmp);
	// cout<<comparator->decrypt_index(ctxt_vec[0])<<endl;
	comparator->print_decrypted(ctxt_vec[0]);
	comparator->print_decrypted(ctxt_vec[2]);

	ctxt_vec[0]-=ctxt_vec[1];
	ctxt_vec[2]-=ctxt_vec[1];
	comparator->compare(ctxt_vec[4],ctxt_vec[0],ctxt_vec[2]);
	cout<<comparator->dec_compare_res(ctxt_vec[4])<<endl;
	comparator->print_decrypted(ctxt_vec[0]);
	comparator->print_decrypted(ctxt_vec[2]);

	// cout<<comparator->decrypt_index(ctxt_vec[0]);

	// clock_t start, end;
	// start = clock();
	// // ctxt_vec[0].multiplyBy(ctxt_vec[1]);
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];
	// ctxt_vec[0] *= ctxt_vec[1];

	// end = clock();
	// cout<<comparator->decrypt_index(ctxt_vec[0])<<endl;

	// cout<<std::fixed<<(double)(end - start)/CLOCKS_PER_SEC<<"s"<<endl;
	// comparator->compare(ctxt_vec[2], ctxt_vec[1], ctxt_vec[0]);
	// ctxt_vec[1] -= ctxt_vec[2];
	// comparator->print_decrypted(ctxt_vec[3]);
	// ctxt_vec[3] *= 1l;
	// cout<<"mul 1"<<endl;
	// comparator->print_decrypted(ctxt_vec[3]);

	// comparator->compare(ctxt_vec[2], ctxt_vec[3], ctxt_vec[1]);
	// comparator->print_decrypted(ctxt_vec[2]);

	// vector<Ctxt>result = comparator->min_dist(ctxt_vec, ctxt_one);
	// // cout<<comparator.decrypt_index(index);


	// for(int i=0; i<result.size(); i++)
	// 	cout<<comparator->decrypt_index(result[i])<<endl;
	// 负数的比较是可行的
	return 0;

}