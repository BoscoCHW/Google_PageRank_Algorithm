//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include <iostream>
#include <vector>
#include "matrix.hpp"
using namespace std;
int main() {
    vector<double> v1{1,2,3,4};
    Matrix m{v1};
    Matrix b{4,2};
    b.setValue(0,0,1);
    b.setValue(0,1,2);
    b.setValue(1,0,3);

    m *= b;
    cout << m << endl;
    cout << b << endl;
    return 0;
}
