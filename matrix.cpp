//
// Created by bosco on 21/09/22.
//

#include "matrix.hpp"

Matrix::Matrix(int n) {
    if (n < 1) {
        throw invalid_argument("Please provide a positive integer.");
    }
    matrix = vector(n, vector<double>(n, 0.0));
}

Matrix::Matrix(int r, int c) {
    if (r < 1 || c < 1) {
        throw invalid_argument("Please provide positive integers.");
    }
    matrix = vector(r, vector<double>(c, 0.0));
}

Matrix::Matrix(const vector<double> &v) {
    unsigned long vectorSize = v.size();
    double squareRoot = sqrt(vectorSize);

    if (squareRoot - floor(squareRoot) == 0) {
        matrix = vector(vectorSize, v);
    } else {
        throw invalid_argument("Please provide a vector with a size of an integer square root.");
    }
}

Matrix::Matrix(vector<vector<double>> m): matrix(std::move(m)){
}

Matrix::Matrix(const Matrix& m) {
    MatrixSize matrixSize = m.getSize();

    this->matrix = vector(matrixSize.rows, vector<double>(matrixSize.cols));
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] = m.matrix.at(i).at(j); // throw out of range error somehow!
        }
    }
//    copy(m.matrix.begin(), m.matrix.end(), back_inserter(matrix));

//    for (const vector<double> &v : m.matrix) {
//        this->matrix.push_back(v);
//    }

}

Matrix::~Matrix() = default;

MatrixSize Matrix::getSize() const {
    int rows = matrix.size();
    int cols = matrix.at(0).size();
    MatrixSize ms{rows, cols};
    return ms;
}

void Matrix::setValue(int row, int col, double val) {
    MatrixSize matrixSiz = getSize();
    if (row < 0 || row >= matrixSiz.rows || col < 0 || col >= matrixSiz.cols) {
        throw invalid_argument("Index out of bound.");
    }
    matrix[row][col] = val;
}
double Matrix::getValue(int row, int col) const {
    MatrixSize matrixSiz = getSize();
    if (row < 0 || row >= matrixSiz.rows || col < 0 || col >= matrixSiz.cols) {
        throw invalid_argument("Index out of bound.");
    }
    return matrix.at(row).at(col);
}

void Matrix::clear() {
    MatrixSize matrixSiz = getSize();
    for (int i = 0; i < matrixSiz.rows; ++i) {
        for (int j = 0; j < matrixSiz.cols; ++j) {
            matrix[i][j] = 0.0;
        }
    }
}

ostream & operator<< (ostream &out, const Matrix &m) {
    MatrixSize matrixSize = m.getSize();
    out << "[";
    for (int i = 0; i < matrixSize.rows; ++i) {
        vector v = m.matrix.at(i);
        out << "[";
        for (int j = 0; j < matrixSize.cols; ++j) {
            out << v.at(j);
            if (j != v.size() - 1) {
                out << ", ";
            }
        }
        if (i != matrixSize.rows - 1) {
            out << "],\n";
        } else {
            out << "]";
        }
    }
    out << "]\n";
    return out;
}

bool operator==(const Matrix &lhs, const Matrix &rhs) {
    MatrixSize matrixSize1 = lhs.getSize();
    MatrixSize matrixSize2 = rhs.getSize();
    if (matrixSize1.rows != matrixSize2.rows || matrixSize1.cols != matrixSize2.cols) {
        return false;
    }
    for (int i = 0; i < matrixSize1.rows; ++i) {
        for (int j = 0; j < matrixSize2.cols; ++j) {
            if (std::abs(lhs.matrix.at(i).at(j) - lhs.matrix.at(i).at(j)) < TOLERANCE) {
                return false;
            }
        }
    }
    return true;
}

Matrix &Matrix::operator++() {
    MatrixSize matrixSize = this->getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] += 1;
        }
    }
    return *this;
}

Matrix Matrix::operator++(int) {
    Matrix temp(*this);
    operator++();
    return temp;
}

Matrix &Matrix::operator--() {
    MatrixSize matrixSize = this->getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] -= 1;
        }
    }
    return *this;
}

Matrix Matrix::operator--(int) {
    Matrix temp(*this);
    operator--();
    return temp;
}

// somehow not called on assignment!
Matrix &Matrix::operator=(Matrix rhs) {
    mySwap(*this, rhs);
    return *this;
}

void mySwap(Matrix &m1, Matrix &m2) {
    m1.matrix.clear();
//    cout << m1 << endl;
    MatrixSize matrixSize = m2.getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        m1.matrix.push_back(m2.matrix.at(i));
    }
}

Matrix &Matrix::operator*=(const Matrix &rhs) {
    MatrixSize mySize = this->getSize();
    MatrixSize rhsSize = rhs.getSize();
    if (mySize.cols != rhsSize.rows) {
        throw invalid_argument("Number of rows of rhs should match number of cols of lhs");
    }

    vector<vector<double>> m;
    for (int i = 0; i < mySize.rows; ++i) {
        vector<double> v;
        for (int j = 0; j < rhsSize.cols; ++j) {
            double sumOfProd = 0;
            for (int p = 0; p < rhsSize.rows; ++p) {
                sumOfProd += rhs.matrix[p][j] * this->matrix[i][p];
            }
            v.push_back(sumOfProd);
        }
        m.push_back(v);
    }

    this->matrix = m;
    return *this;
}

Matrix operator*(Matrix lhs, const Matrix &rhs) {
    lhs *= rhs;
    return lhs;
}

Matrix &Matrix::operator+=(const Matrix &rhs) {
    MatrixSize ms1 = this->getSize();
    MatrixSize ms2 = rhs.getSize();
    if (ms1.rows != ms2.rows || ms1.cols != ms2.cols) {
        throw invalid_argument("Different sizes of matrix.");
    }
    for (int i = 0; i < ms1.rows; ++i) {
        for (int j = 0; j < ms1.cols; ++j) {
            matrix[i][j] += rhs.matrix[i][j];
        }
    }
    return *this;
}

Matrix operator+(Matrix lhs, const Matrix &rhs) {
    lhs += rhs;
    return lhs;
}

Matrix &Matrix::operator-=(const Matrix &rhs) {
    MatrixSize ms1 = this->getSize();
    MatrixSize ms2 = rhs.getSize();
    if (ms1.rows != ms2.rows || ms1.cols != ms2.cols) {
        throw invalid_argument("Different sizes of matrix.");
    }
    for (int i = 0; i < ms1.rows; ++i) {
        for (int j = 0; j < ms1.cols; ++j) {
            matrix[i][j] -= rhs.matrix[i][j];
        }
    }
    return *this;
}

Matrix operator-(Matrix lhs, const Matrix &rhs) {
    lhs -= rhs;
    return lhs;
}

Matrix &Matrix::operator*=(const double rhs) {
    MatrixSize ms = getSize();
    for (int i = 0; i < ms.rows; ++i) {
        for (int j = 0; j < ms.cols; ++j) {
            matrix[i][j] *= rhs;
        }
    }
    return *this;
}

Matrix operator*(Matrix lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
}
