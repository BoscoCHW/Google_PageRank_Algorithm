//
// Created by bosco on 21/09/22.
//

#include "matrix.hpp"

Matrix::Matrix(int n) {
    if (n < 1) {
        throw invalid_argument("Please provide a positive integer.");
    }
    matrix = vector(n, vector<double>(n));
}

Matrix::Matrix(int r, int c) {
    if (r < 1 || c < 1) {
        throw invalid_argument("Please provide positive integers.");
    }
    matrix = vector(r, vector<double>(c));
}

Matrix::Matrix(vector<double> v) {
    int vectorSize = v.size();
    double squareRoot = sqrt(vectorSize);

    if (squareRoot - floor(squareRoot) == 0) {
        matrix = vector(vectorSize, v);
    } else {
        throw invalid_argument("Please provide a vector with a size of an integer square root.");
    }
}

Matrix::Matrix(const Matrix& m) {
    copy(m.matrix.begin(), m.matrix.end(), back_inserter(matrix));
}

tuple<int, int> Matrix::getSize() {
    int rows = matrix.size();
    int cols = matrix.at(0).size();
    return {rows, cols}
}

void Matrix::setValue(int row, int col, double val) {
    auto [rows, cols] = getSize();
    if (row >= rows || col >= cols) {
        throw invalid_argument("Index out of bound.");
    }
    matrix.at(row).at(col) = val;
}
double Matrix::getValue(int row, int col) {
    auto [rows, cols] = getSize();
    if (row >= rows || col >= cols) {
        throw invalid_argument("Index out of bound.");
    }
    return matrix.at(row).at(col);
}

void Matrix::clear() {
    matrix.
}