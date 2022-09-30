//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include <iostream>
#include <vector>
#include "matrix.hpp"
using namespace std;
int main() {
    vector<double> v{2.2, 3, 4, 5, 2, 3, 4, 5, 9};
    Matrix m{v};
    std::cout << m << std::endl;
    m.clear();
    cout << m << endl;
    m.setValue(2,4 ,100);
    cout << m << endl;
    return 0;
}
