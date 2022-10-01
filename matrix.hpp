//
// Created by bosco on 21/09/22.
//

#ifndef LAB1TEMPLATE_MATRIX_HPP
#define LAB1TEMPLATE_MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <tuple>
#include <cmath>
#include <iostream>

using namespace std;

constexpr double TOLERANCE = 0.0001;

struct MatrixSize {
    int rows;
    int cols;
};

class Matrix {
private:
    vector<vector<double>> matrix;
public:
    Matrix(): matrix{{0}} {};
    Matrix(int n);
    Matrix(int r, int c);
    Matrix(vector<double> v);
    Matrix(const Matrix& m);
    ~Matrix();
    MatrixSize getSize() const;
    void setValue(int row, int col, double val);
    double getValue(int row, int col) const;
    void clear();
    friend ostream & operator << (ostream &out, const Matrix &m);
    friend bool operator== (const Matrix& lhs, const Matrix& rhs);
    friend bool operator!= (const Matrix& lhs, const Matrix& rhs) { return !(lhs == rhs); }
    Matrix& operator++();
    Matrix operator++(int);
    Matrix& operator--();
    Matrix operator--(int);
    Matrix& operator=(Matrix rhs);
    void mySwap(Matrix& m1, Matrix& m2);
};


#endif //LAB1TEMPLATE_MATRIX_HPP
