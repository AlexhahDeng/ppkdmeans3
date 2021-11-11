#include <iostream>

#include "my_tools.h"

// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
// B 7 2 300 170 3 10 y
int main(){
	Comparator *comparator = generate_comparator(false);

	int runs = 1;
	vector<long int>tmp{5,4,8,2,5};
	Ctxt ctxt_one = comparator->gen_ctxt_one();

	vector<Ctxt>ctxt_vec = comparator->encrypt_vector(tmp);
	// comparator->compare(ctxt_vec[2], ctxt_vec[1], ctxt_vec[0]);

	vector<Ctxt>result = comparator->min_dist(ctxt_vec, ctxt_one);
	// cout<<comparator.decrypt_index(index);


	for(int i=0; i<result.size(); i++)
		cout<<comparator->decrypt_index(result[i])<<endl;


	return 0;

}