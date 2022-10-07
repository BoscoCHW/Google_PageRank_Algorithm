//
// Created by Bosco on 06-Oct.-2022.
//

#include "pageRanker.hpp"

/**
 * PageRanker constructor
 * @param connectivityMatrix
 */
PageRanker::PageRanker(const Matrix& connectivityMatrix) {
    MatrixSize matrixSize = connectivityMatrix.getSize();
    if (matrixSize.rows != matrixSize.cols) {
        throw invalid_argument("Connectivity matrix must be a square matrix");
    }
    this->connectivityMatrix = connectivityMatrix;
}

/**
 * Create a connectivity matrix from a text file
 * @param filename
 */
void PageRanker::loadConnectivityMatrixFromFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }
    string line;
    int num;
    vector<double> nums;
    while (getline(file, line)) {
        istringstream iss(line);
        while (iss >> num) {
            nums.push_back(num);
        }
    }
    Matrix matrix(nums);
    connectivityMatrix = matrix;
}

/**
 * Pre-condition: Connectivity matrix loaded
 * Rank pages using the Google PageRank algorithm
 * @return rank
 */
Matrix PageRanker::rankPages() {
    Matrix stochasticMatrix = createStochasticMatrix();
    Matrix probabilityMatrix = createProbabilityMatrix(stochasticMatrix);
    MatrixSize matrixSize = probabilityMatrix.getSize();

    Matrix rank(matrixSize.rows, 1);
    for (int i = 0; i < matrixSize.rows; ++i) {
        rank.setValue(i, 0, 1.0);
    }

    while (true) {
        Matrix newRank = probabilityMatrix * rank;
        if (abs(newRank.getValue(0, 0) - rank.getValue(0, 0)) < PageRanker::MARKOV_THRESHOLD) {
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
 * From connectivity matrix, create a stochastic matrix with normalized probability
 * @return stochastic matrix
 */
Matrix PageRanker::createStochasticMatrix() {
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
Matrix PageRanker::createProbabilityMatrix(Matrix &stochasticMatrix) {
    MatrixSize matrixSize = stochasticMatrix.getSize();

    vector<double> evenDistribution(matrixSize.rows * matrixSize.cols);
    for (int i = 0; i < matrixSize.rows * matrixSize.cols; ++i) {
        evenDistribution[i] = 1.0 / matrixSize.rows;
    }

    constexpr double p = 0.85;

    Matrix Q(evenDistribution);

    return stochasticMatrix * p + Q * (1 - p);
}
