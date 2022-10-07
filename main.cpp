//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "matrix.hpp"
using namespace std;

/**
 * Threshold used in the Markov Process.
 */
constexpr double MARKOV_THRESHOLD = 0.0001;

/**
 * Create a connectivity matrix from a text file
 * @param filename
 * @return a connectivity matrix
 */
Matrix createConnectivityMatrixFromFile(const string filename) {
    ifstream connectivity(filename);
    if (!connectivity.is_open()) {
        cout << "Cannot open file." << endl;
    }
    string line;
    int num;
    vector<double> nums;
    while (getline(connectivity, line)) {
        istringstream iss(line);
        while (iss >> num) {
            nums.push_back(num);
        }
    }
    Matrix matrix(nums);

    return matrix;
}

/**
 * Given a connectivity matrix, create a stochastic matrix
 * with normalized probability
 * @param connectivityMatrix
 * @return stochastic matrix
 */
Matrix createStochasticMatrix(const Matrix &connectivityMatrix) {
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
    return stochasticMatrix;
}

/**
 * Create a probability matrix from stochastic matrix
 * Factor in random walk probability
 * @param stochasticMatrix
 * @return probability matrix
 */
Matrix createProbabilityMatrix(Matrix &stochasticMatrix) {
    MatrixSize matrixSize = stochasticMatrix.getSize();

    vector<double> evenDistribution;
    for (int i = 0; i < matrixSize.rows * matrixSize.cols; ++i) {
        evenDistribution.push_back(1.0 / matrixSize.rows);
    }

    constexpr double p = 0.85;

    Matrix Q(evenDistribution);

    return stochasticMatrix * p + Q * (1 - p);
}

/**
 * Given the probability matrix for a web,
 * rank pages using the Google PageRank algorithm
 * @param probabilityMatrix
 * @return rank
 */
Matrix rankPages(Matrix &probabilityMatrix) {
    MatrixSize matrixSize = probabilityMatrix.getSize();

    Matrix rank(matrixSize.rows, 1);
    for (int i = 0; i < matrixSize.rows; ++i) {
        rank.setValue(i, 0, 1.0);
    }

    while (true) {
        Matrix newRank = probabilityMatrix * rank;
        if (abs(newRank.getValue(0, 0) - rank.getValue(0, 0)) < MARKOV_THRESHOLD) {
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
    return rank;
}

/**
 * Read a file named "connectivity.txt" containing the connectivity matrix of a web
 * Rank each page using the Google PageRank algorithm.
 */
int main() {

    Matrix connectivityMatrix = createConnectivityMatrixFromFile("../connectivity.txt");

    Matrix stochasticMatrix = createStochasticMatrix(connectivityMatrix);

    Matrix probabilityMatrix = createProbabilityMatrix(stochasticMatrix);

    Matrix rank = rankPages(probabilityMatrix);

    MatrixSize rankSize = rank.getSize();
    char pageName = 'A';
    for (int i = 0; i < rankSize.rows; ++i) {
        double pageRank = rank.getValue(i, 0) * 100;
        printf("Page %c: %0.2f%%\n", pageName++, pageRank);
    }

    return 0;
}
