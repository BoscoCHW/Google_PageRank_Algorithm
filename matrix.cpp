//
// Created by bosco on 21/09/22.
//

#include "matrix.hpp"

/**
 * Matrix constructor:
 * Create a matrix of size n * n
 * @param n
 */
Matrix::Matrix(const int n) {
    if (n < 1) {
        throw invalid_argument("Please provide a positive integer.");
    }
    matrix = vector(n, vector<double>(n, 0.0));
}

/**
 * Matrix constructor:
 * Create a matrix of size r * c
 * @param r number of rows
 * @param c number of columns
 */
Matrix::Matrix(const int r, const int c) {
    if (r < 1 || c < 1) {
        throw invalid_argument("Please provide positive integers.");
    }
    matrix = vector(r, vector<double>(c, 0.0));
}

/**
 * Matrix constructor:
 * Create a matrix from a vector, reshaped into a square matrix
 * @param v vector of size of an square number
 */
Matrix::Matrix(const vector<double> &v) {
    unsigned long vectorSize = v.size();
    double squareRoot = sqrt(vectorSize);

    if (squareRoot - floor(squareRoot) != 0) {
        throw invalid_argument("Please provide a vector with a size of an integer square root.");
    }
    int matrixDimension = (int) squareRoot;
    vector<vector<double>> m(matrixDimension);
    for (int i = 0; i < matrixDimension; ++i) {
        vector<double> vector(matrixDimension);
        for (int j = 0; j < matrixDimension; ++j) {
            vector[j] = v.at(i * matrixDimension + j);
        }
        m[i] = vector;
    }

    matrix = m;
}

/**
 * Matrix copy constructor
 * @param m
 */
Matrix::Matrix(const Matrix& m) {
    MatrixSize matrixSize = m.getSize();

    this->matrix = vector(matrixSize.rows, vector<double>(matrixSize.cols));
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] = m.matrix.at(i).at(j);
        }
    }

}

/**
 * Matrix destructor
 */
Matrix::~Matrix() = default;

/**
 * Return the size of the matrix
 * @return MatrixSize
 */
MatrixSize Matrix::getSize() const {
    int rows = matrix.size();
    int cols = matrix.at(0).size();
    MatrixSize ms{rows, cols};
    return ms;
}

/**
 * Set value to a specified cell of the matrix
 * @param row
 * @param col
 * @param val
 */
void Matrix::setValue(int row, int col, double val) {
    MatrixSize matrixSiz = getSize();
    if (row < 0 || row >= matrixSiz.rows || col < 0 || col >= matrixSiz.cols) {
        throw invalid_argument("Index out of bound.");
    }
    matrix[row][col] = val;
}

/**
 * Get value of a cell in the matrix
 * @param row
 * @param col
 * @return the value at matrix[row][col]
 */
double Matrix::getValue(int row, int col) const {
    MatrixSize matrixSiz = getSize();
    if (row < 0 || row >= matrixSiz.rows || col < 0 || col >= matrixSiz.cols) {
        throw invalid_argument("Index out of bound.");
    }
    return matrix.at(row).at(col);
}

/**
 * Change value of all cells to 0.0
 */
void Matrix::clear() {
    MatrixSize matrixSiz = getSize();
    for (int i = 0; i < matrixSiz.rows; ++i) {
        for (int j = 0; j < matrixSiz.cols; ++j) {
            matrix[i][j] = 0.0;
        }
    }
}

/**
 * Overload insertion operator
 * Print information of the matrix
 * @param out
 * @param m
 * @return ostream
 */
ostream & operator<< (ostream &out, const Matrix &m) {
    MatrixSize matrixSize = m.getSize();
    out << "[";
    for (int i = 0; i < matrixSize.rows; ++i) {
        vector v = m.matrix.at(i);
        out << "[";
        for (int j = 0; j < matrixSize.cols; ++j) {
            out << v.at(j);
            if (j != (int)v.size() - 1) {
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

/**
 * Overload equals operator
 * @param lhs
 * @param rhs
 * @return true if both matrices have the same size and same values in each cell
 */
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

/**
 * Overload prefix increment operator
 * Increase value of each cell by 1
 */
Matrix &Matrix::operator++() {
    MatrixSize matrixSize = this->getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] += 1;
        }
    }
    return *this;
}

/**
 * Overload postfix increment operator
 */
Matrix Matrix::operator++(int) {
    Matrix temp(*this);
    operator++();
    return temp;
}

/**
 * Overload prefix decrement operator
 * Decrease value of each cell by 1
 */
Matrix &Matrix::operator--() {
    MatrixSize matrixSize = this->getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            this->matrix[i][j] -= 1;
        }
    }
    return *this;
}

/**
 * Overload postfix decrement operator
 */
Matrix Matrix::operator--(int) {
    Matrix temp(*this);
    operator--();
    return temp;
}

/**
 * Overload assignment constructor
 */
Matrix &Matrix::operator=(Matrix rhs) {
    mySwap(*this, rhs);
    return *this;
}

/**
 * Assign value of each cell of matrix 2 to matrix 1
 * @param m1 Matrix 1
 * @param m2 Matrix 2
 */
void mySwap(Matrix &m1, Matrix &m2) {
    m1.matrix.clear();
    MatrixSize matrixSize = m2.getSize();
    for (int i = 0; i < matrixSize.rows; ++i) {
        m1.matrix.push_back(m2.matrix.at(i));
    }
}

/**
 * Overload *= operator
 * Perform matrix multiplication
 */
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

/**
 * Overload * operator
 * Perform matrix multiplication
 */
Matrix operator*(Matrix lhs, const Matrix &rhs) {
    lhs *= rhs;
    return lhs;
}

/**
 * Overload += operator
 * Perform matrix addition
 */
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

/**
 * Overload + operator
 */
Matrix operator+(Matrix lhs, const Matrix &rhs) {
    lhs += rhs;
    return lhs;
}

/**
 * Overload -= operator
 * Perform matrix subtraction
 */
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

/**
 * Overload - operator
 */
Matrix operator-(Matrix lhs, const Matrix &rhs) {
    lhs -= rhs;
    return lhs;
}

/**
 * Overload *= operator for scalar multiplication
 */
Matrix &Matrix::operator*=(const double rhs) {
    MatrixSize ms = getSize();
    for (int i = 0; i < ms.rows; ++i) {
        for (int j = 0; j < ms.cols; ++j) {
            matrix[i][j] *= rhs;
        }
    }
    return *this;
}

/**
 * Overload * operator for scalar multiplication
 */
Matrix operator*(Matrix lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
}
