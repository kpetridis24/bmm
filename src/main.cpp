/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <triCounting.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main()
{
    /* ------------------------------- read matrix ------------------------------ */

    int n;
    int nnz;

    std::string graph = "s6.mtx";
    std::string file = "graphs/" + graph;

    readMtxValues(file, n, nnz);

    coo M;
    util::initCoo(M, n, nnz);

    openMtxFile(file, M.col, M.row, M.n, M.nnz);

    csr A;
    util::initCsr(A, n, nnz);
    csc B;
    util::initCsc(B, n, nnz);

    // std::cout << "\nCOO row:\t";
    // prt::arr(M.row, nnz);
    // std::cout << "COO col:\t";
    // prt::arr(M.col, nnz);

    coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
    coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);

    util::delCoo(M);

    // std::cout << "\nCSR row_ptr:";
    // prt::arr(A.rowPtr, n + 1);
    // std::cout << "CSR col_ind:";
    // prt::arr(A.colInd, nnz);

    // std::cout << "\nCSC col_ptr:";
    // prt::arr(B.colPtr, n + 1);
    // std::cout << "CSC row_ind:";
    // prt::arr(B.rowInd, nnz);

    /* -------------------------------- bmm test -------------------------------- */

    coo C;
    util::initCoo(C, A.n, A.nnz * B.nnz); // TODO check max size
    // bmm(A, B, C);
    maskedBmm(A, A, B, C);

    prt::arr(C.row, C.nnz);
    prt::arr(C.col, C.nnz);

    /* ------------------------------- free memory ------------------------------ */

    util::delCsr(A);
    util::delCsc(B);
    util::delCoo(C);

    return 0;

}