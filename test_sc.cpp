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

	vector<long int>tmp{100,10,9,255,5,117645};

	vector<Ctxt>ctxt_vec = comparator->encrypt_vector(tmp);
	vector<int>ss1(ctxt_vec.size()), ss2(ss1);

	// ctxt_vec[0] -= ctxt_vec[1];
	// ctxt_vec[2] -= ctxt_vec[1];
	// comparator->print_decrypted(ctxt_vec[0]);
	// comparator->print_decrypted(ctxt_vec[2]);
	
	// comparator->compare(lessthan, ctxt_vec[0], ctxt_vec[2]);
	// comparator->print_decrypted(lessthan);

	comparator->he_to_ss(ctxt_vec, ss1, ss2);
	

	return 0;

}