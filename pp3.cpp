#include <iostream>
#include "cloud.h"
#include "func.h"
#include "tools.h"

using namespace std;

int main(){
    Comparator* comparator = generate_comparator();
    comparator->test_compare();
    return 0;
}