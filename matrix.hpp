//
// Created by bosco on 21/09/22.
//

#ifndef LAB1TEMPLATE_MATRIX_HPP
#define LAB1TEMPLATE_MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>

using namespace std;

constexpr double TOLERANCE = 0.0001;

/**
 * A container for the size of a matrix
 */
struct MatrixSize {
    int rows;
    int cols;
};

class Matrix {
private:
    vector<vector<double>> matrix;
public:
    Matrix(): matrix{{0}} {};
    explicit Matrix(int n);
    Matrix(int r, int c);
    explicit Matrix(const vector<double> &v);
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
    friend void mySwap(Matrix& m1, Matrix& m2);

    Matrix& operator+=(const Matrix& rhs);
    friend Matrix operator+(Matrix lhs, const Matrix& rhs);

    Matrix& operator-=(const Matrix& rhs);
    friend Matrix operator-(Matrix lhs, const Matrix& rhs);

    Matrix& operator*=(const Matrix& rhs);
    friend Matrix operator*(Matrix lhs, const Matrix& rhs);

    Matrix& operator*=(const double rhs);
    friend Matrix operator*(Matrix lhs, const double rhs);
};


#endif //LAB1TEMPLATE_MATRIX_HPP
