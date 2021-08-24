/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <unistd.h>

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <block-bmm.cpp>
#include <triCounting.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main()
{
    struct timeval timer;

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

    // prt::cooMat(M);

    coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
    coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);

    util::delCoo(M);

    // prt::csrMat(A);
    // prt::cscMat(B);

    std::cout << "\nMatrix read successfully\nn = " << A.n << ", nnz = " << A.nnz << std::endl;
    
    /* --------------------------- bcsr blocking test --------------------------- */

    std::cout << "\nBlocking A in B-CSR...\n";

    int b = 2;
    int numBlocks = (n / b) * (n / b);
    int LL_bRowPtrSize = numBlocks * (b + 1);

    bcsr blA;
    blA.n = A.n;
    blA.b = b;

    // init Low-Level CSR
    blA.LL_bRowPtr = new int[LL_bRowPtrSize]();
    blA.LL_bColInd = new int[nnz]();

    // blocking
    ret _ret = csr2bcsr(A, blA);

    blA.HL_bRowPtr = _ret.ret1;
    blA.HL_bColInd = _ret.ret2;
    blA.nzBlockIndex = _ret.ret3;
    blA.blockNnzCounter = _ret.ret4;

    /* --------------------------- bcsc blocking test --------------------------- */

    std::cout << "\nBlocking B in B-CSC...\n";

    b = 2;
    int LL_bColPtrSize = numBlocks * (b + 1);

    bcsc blB;
    blB.n = A.n;
    blB.b = b;

    // init Low-Level CSC
    blB.LL_bColPtr = new int[LL_bColPtrSize]();
    blB.LL_bRowInd = new int[nnz]();

    // blocking
    _ret = csr2bcsc(B, blB);

    blB.HL_bColPtr = _ret.ret1;
    blB.HL_bRowInd = _ret.ret2;
    blB.nzBlockIndex = _ret.ret3;
    blB.blockNnzCounter = _ret.ret4;

    /* -------------------------------- bmm test -------------------------------- */

    // coo C;

    // timer = util::tic();
    
    // // util::initCoo(C, A.n, A.nnz * B.nnz); // TODO check max size
    // // bmm(A, B, C);

    // util::initCoo(C, A.n, A.nnz);
    // maskedBmm(A, A, B, C);

    // // util::initCoo(C, A.n, A.nnz);
    // // maskedBmm2(A, A, A, C);

    // // prt::cooMat(C);

    // double t = util::toc(timer);
    // std::cout << "BMM completed\nC.nnz = " << C.nnz << "\nBMM time = " << t << " seconds" << std::endl;

    /* ------------------------------- free memory ------------------------------ */

    util::delCsr(A);
    util::delCsc(B);
    // util::delCoo(C);
    util::delBcsr(blA);
    util::delBcsc(blB);

    return 0;
}