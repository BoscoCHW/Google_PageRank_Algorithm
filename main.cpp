//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "matrix.hpp"
using namespace std;

Matrix createConnectivityMatrix(const string filename) {
    ifstream connectivity(filename);
    if (!connectivity.is_open()) {
        cout << "Cannot open file." << endl;
    }
    string line;
    int num;
    vector<vector<double>> m;
    while (getline(connectivity, line)) {
        istringstream iss(line);
        vector<double> v;
        while (iss >> num) {
            v.push_back(num);
        }
        m.push_back(v);
    }
    Matrix matrix(m);

    return matrix;
}

int main() {

    Matrix connectivityMatrix = createConnectivityMatrix("connectivityMatrix.txt");
    MatrixSize matrixSize = connectivityMatrix.getSize();
    Matrix stochasticMatrix(matrixSize.rows, matrixSize.cols);

    vector<double> outDegrees;
    for (int j = 0; j < matrixSize.cols; ++j) {
        double columnSum = 0;
        for (int i = 0; i < matrixSize.rows; ++i) {
            columnSum += connectivityMatrix.getValue(i, j);
        }
        outDegrees.push_back(columnSum);
    }

    for (int i = 0; i < matrixSize.rows; ++i) {
        for (int j = 0; j < matrixSize.cols; ++j) {
            if (std::abs(outDegrees.at(j) - 0) > 0.00001) {
                // if out degree is not 0
                double num = connectivityMatrix.getValue(i, j) / outDegrees.at(j);
                stochasticMatrix.setValue(i, j, num);
            } else {
                // if out degree is 0
                stochasticMatrix.setValue(i, j, 1.0 / matrixSize.rows);
            }

        }
    }
    vector<double> evenDistribution;
    for (int i = 0; i < matrixSize.rows; ++i) {
        evenDistribution.push_back(1.0 / matrixSize.rows);
    }

    double p = 0.85;

    Matrix Q(evenDistribution);
    Matrix probabilityMatrix = stochasticMatrix * p + Q * (1 - p);
    cout << probabilityMatrix << endl;

    Matrix rank(matrixSize.rows, 1);
    for (int i = 0; i < matrixSize.rows; ++i) {
        rank.setValue(i, 0, 1.0);
    }

    double threshold = 0.001;

    while (true) {
        Matrix newRank = probabilityMatrix * rank;
        if (abs(newRank.getValue(0, 0) - rank.getValue(0, 0)) < threshold) {
            rank = newRank;
            break;
        }
        rank = newRank;
    }

    double sumOfRank = 0;

    for (int i = 0; i < matrixSize.rows; ++i) {
        sumOfRank += rank.getValue(i, 0);
    }

    for (int i = 0; i < matrixSize.rows; ++i) {
        rank.setValue(i, 0, rank.getValue(i, 0) / sumOfRank);
    }
    cout << rank << endl;
    return 0;
}
