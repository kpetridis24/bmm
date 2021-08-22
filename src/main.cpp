/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <vector>

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <triCounting.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main()
{
    struct timeval timer;

    /* ------------------------------- read matrix ------------------------------ */

    int n;
    int nnz;

    std::string graph = "com-Youtube.mtx";
    std::string file = "graphs/" + graph;

    readMtxValues(file, n, nnz);

    coo M;
    util::initCoo(M, n, nnz);

    openMtxFile(file, M.col, M.row, M.n, M.nnz);

    csr A;
    util::initCsr(A, n, nnz);
    csc B;
    util::initCsc(B, n, nnz);

    // prt::cooMat(M);

    coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
    coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);

    util::delCoo(M);

    // prt::csrMat(A);
    // prt::cscMat(B);

    std::cout << "\nMatrix read successfully\nn = " << A.n << ", nnz = " << A.nnz << std::endl;
    
    /* -------------------------------- bmm test -------------------------------- */

    coo C;

    timer = util::tic();

    // util::initCoo(C, A.n, A.nnz * B.nnz); // TODO check max size
    // bmm(A, B, C);
    util::initCoo(C, A.n, A.nnz);
    maskedBmm(A, A, B, C);

    // prt::cooMat(C);

    double t = util::toc(timer);
    std::cout << "BMM completed\nC.nnz = " << C.nnz << "\nBMM time = " << t << " seconds" << std::endl;

    /* ------------------------------- free memory ------------------------------ */

    util::delCsr(A);
    util::delCsc(B);
    util::delCoo(C);

    return 0;

}