//
// Created by bosco on 21/09/22.
//

#ifndef LAB1TEMPLATE_MATRIX_HPP
#define LAB1TEMPLATE_MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <tuple>
#include <cmath>

using namespace std;

class Matrix {
private:
    vector<vector<double>> matrix;
public:
    Matrix(): matrix{{0}} {};
    Matrix(int n);
    Matrix(int r, int c);
    Matrix(vector<double> v);
    Matrix(const Matrix& m);
    tuple<int, int> getSize();
    void setValue(int row, int col, double val);
    double getValue(int row, int col);
    void clear();
};


#endif //LAB1TEMPLATE_MATRIX_HPP
