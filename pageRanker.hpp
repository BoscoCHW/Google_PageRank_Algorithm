//
// Created by Bosco on 06-Oct.-2022.
//

#ifndef LAB1TEMPLATE_PAGERANKER_HPP
#define LAB1TEMPLATE_PAGERANKER_HPP

#include <fstream>
#include <sstream>
#include "matrix.hpp"

class PageRanker {
private:
    Matrix connectivityMatrix;
    Matrix createStochasticMatrix();
    Matrix createProbabilityMatrix(Matrix &stochasticMatrix);

public:
    /**
     * Threshold used in the Markov Process.
     */
    constexpr static double MARKOV_THRESHOLD = 0.0001;
    PageRanker() = default;
    explicit PageRanker(const Matrix& connectivityMatrix);
    void loadConnectivityMatrixFromFile(const string& filename);
    Matrix rankPages();

};


#endif //LAB1TEMPLATE_PAGERANKER_HPP
