//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include <iostream>
#include <vector>
#include "matrix.hpp"
using namespace std;
int main() {
    vector<double> v{2.2, 3, 4, 5, 2, 3, 4, 5, 9};
    Matrix m{v};

    Matrix b{3};
    b.setValue(0,0,9);
    b.setValue(0,1,8);
    m = b;
    cout << m << endl;

    m.setValue(0,1,100);
    cout << m << endl;
    cout << b << endl;
    return 0;
}
