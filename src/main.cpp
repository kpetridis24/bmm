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

    int *cooRow = new int[nnz]();
    int *cooCol = new int[nnz]();

    openMtxFile(file, cooCol, cooRow, n, nnz);

    csr A;
    util::initCsr(&A, n, nnz);
    csc B;
    util::initCsc(&B, n, nnz);

    // std::cout << "\nCOO row:\t";
    // prt::arr(coo_row, nnz);
    // std::cout << "COO col:\t";
    // prt::arr(coo_col, nnz);

    coo2csr(A.rowPtr, A.colInd, cooRow, cooCol, A.nnz, A.n, 0);
    coo2csr(B.colPtr, B.rowInd, cooCol, cooRow, B.nnz, B.n, 0);

    delete[] cooRow;
    delete[] cooCol;

    std::cout << "\nCSR row_ptr:";
    prt::arr(A.rowPtr, n + 1);
    std::cout << "CSR col_ind:";
    prt::arr(A.colInd, nnz);

    std::cout << "\nCSC col_ptr:";
    prt::arr(B.colPtr, n + 1);
    std::cout << "CSC row_ind:";
    prt::arr(B.rowInd, nnz);

    /* ------------------------------ blocking test ----------------------------- */

    // int b = 2;
    // int numBlocks = (n / b) * (n / b);
    // int LL_bRowPtrSize = numBlocks * (b + 1);
    // int blocksPerRow = n / b;

    // int *nzBlockIndex;
    // int *blockNnzCounter;

    // // Low-Level CSR
    // int *LL_bRowPtr = new int[LL_bRowPtrSize]();
    // int *LL_bColInd = new int[nnz]();

    // // High-Level B-CSR
    // int *HL_bRowPtr;
    // int *HL_bColInd;

    // // blocking
    // ret _ret = csr2blocks(A.rowPtr, A.colInd, n, nnz, b, LL_bRowPtr, LL_bColInd);

    // HL_bRowPtr = _ret.ret1;
    // HL_bColInd = _ret.ret2;
    // nzBlockIndex = _ret.ret3;
    // blockNnzCounter = _ret.ret4;

    /* ---------------------------- triCounting test ---------------------------- */

    // int trNum = bCsrTriCount( LL_bRowPtr, 
    //                         LL_bColInd, 
    //                         HL_bRowPtr, 
    //                         HL_bColInd,
    //                         nzBlockIndex,
    //                         blockNnzCounter, 
    //                         n, 
    //                         b );

    // std::cout << "Num of triangles: " << trNum << std::endl;

    /* -------------------------------- bmm test -------------------------------- */



    /* ------------------------------- free memory ------------------------------ */

    delete[] A.rowPtr;
    delete[] A.colInd;
    // delete[] LL_bRowPtr;
    // delete[] LL_bColInd;
    // delete[] HL_bRowPtr;
    // delete[] HL_bColInd;
    // delete[] nzBlockIndex;
    // delete[] blockNnzCounter;

    return 0;

}