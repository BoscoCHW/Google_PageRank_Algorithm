//Name: Heung Wai Chan (Bosco)
//Student# : A01258687

#include "matrix.hpp"
#include "pageRanker.hpp"
using namespace std;

/**
 * Read a file named "connectivity.txt" containing the connectivity matrix of a web
 * Rank each page using the Google PageRank algorithm.
 */
int main() {

    PageRanker pageRanker;
    pageRanker.loadConnectivityMatrixFromFile("../connectivity.txt");
    Matrix rank = pageRanker.rankPages();

    MatrixSize rankSize = rank.getSize();
    char pageName = 'A';
    for (int i = 0; i < rankSize.rows; ++i) {
        double pageRank = rank.getValue(i, 0) * 100;
        printf("Page %c: %0.2f%%\n", pageName++, pageRank);
    }

    return 0;
}
