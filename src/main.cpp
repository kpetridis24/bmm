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
#include <bits/stdc++.h>
#include <mpi.h>

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <block-bmm.cpp>
#include <masked-block-bmm.cpp>
// #include <parallel-masked-block-bmm.cpp>
#include <distributed-block-bmm.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main(int argc, char **argv)
{   
    /* -------------------------------------------------------------------------- */
    /*                        OpenMPI BMM distribution test                       */
    /* -------------------------------------------------------------------------- */

    struct timeval timer;
    double t = -1;

    /* --------------------------- initialize OpenMPI --------------------------- */

    int numProcesses, rank;
    int graphInd = 2;
    int blockSizeA;
    int blockSizeB;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    coo cooA;
    coo cooB;
    coo _cooA;

    /* ---------------------- distribute A and broadcast B ---------------------- */

    distributeCooMatrix(numProcesses, rank, cooA, _cooA, graphInd, blockSizeA);
    broadcastCooMatrix(numProcesses, rank, cooB, graphInd, blockSizeB);

    /* ---------------- fix indices of matrix A - remove offsets ---------------- */

    int rowsPerChunk = _cooA.n / numProcesses;
    int chunkStartingRow = rowsPerChunk * rank;

    util::removeCooRowOffsets(_cooA, chunkStartingRow);

    /* ---------------------- convert A to CSR and B to CSC --------------------- */

    csr csrA;
    util::initCsr(csrA, _cooA.m, _cooA.n, _cooA.nnz);

    coo2csr(csrA.rowPtr, csrA.colInd, _cooA.row, _cooA.col, _cooA.nnz, _cooA.m, 0);

    // prt::csrMat(csrA);

    csc cscB;
    util::initCsc(cscB, cooB.m, cooB.n, cooB.nnz);

    coo2csr(cscB.colPtr, cscB.rowInd, cooB.col, cooB.row, cooB.nnz, cooB.n, 0);

    // prt::cscMat(cscB);

    /* -------------------------- blocking of A Matrix -------------------------- */

    timer = util::tic();

    bcsr bcsrA;
    bcsrA.m = _cooA.m;
    bcsrA.n = _cooA.n;
    bcsrA.b = blockSizeA;

    int numBlocksA = (bcsrA.m / bcsrA.b) * (bcsrA.n / bcsrA.b);
    int LL_bRowPtrSize = numBlocksA * (bcsrA.b + 1);

    // init Low-Level CSR
    bcsrA.LL_bRowPtr = new int[LL_bRowPtrSize]();
    bcsrA.LL_bColInd = new int[_cooA.nnz]();

    // blocking
    ret _ret1 = csr2bcsr(csrA, bcsrA);

    bcsrA.HL_bRowPtr = _ret1.ret1;
    bcsrA.HL_bColInd = _ret1.ret2;
    bcsrA.nzBlockIndex = _ret1.ret3;
    bcsrA.blockNnzCounter = _ret1.ret4;
    int HL_bRowPtrSize = _ret1.size1;
    int HL_bColIndSize = _ret1.size2;
    int nzBlockIndexSizeA = _ret1.size3;
    int blockNnzCounterSizeA = _ret1.size4;


    t = util::toc(timer);
    std::cout << "\nBlocking A in B-CSR completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    /* -------------------------- blocking of B Matrix -------------------------- */

    timer = util::tic();

    bcsc bcscB;
    bcscB.m = cooB.m;
    bcscB.n = cooB.n;
    bcscB.b = blockSizeB;

    int numBlocks = (bcscB.m / bcscB.b) * (bcscB.n / bcscB.b);
    int LL_bColPtrSize = numBlocks * (bcscB.b + 1);

    // init Low-Level CSC
    bcscB.LL_bColPtr = new int[LL_bColPtrSize]();
    bcscB.LL_bRowInd = new int[cooB.nnz]();

    // blocking
    ret _ret2 = csc2bcsc(cscB, bcscB);
    
    bcscB.HL_bColPtr = _ret2.ret1;
    bcscB.HL_bRowInd = _ret2.ret2;
    bcscB.nzBlockIndex = _ret2.ret3;
    bcscB.blockNnzCounter = _ret2.ret4;
    int HL_bColPtrSize = _ret2.size1;
    int HL_bRowIndSize = _ret2.size2;
    int nzBlockIndexSizeB = _ret2.size3;
    int blockNnzCounterSizeB = _ret2.size4;


    t = util::toc(timer);
    std::cout << "\nBlocking B in B-CSC completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    /* -------------------------------- block BMM ------------------------------- */

    timer = util::tic();

    std::multimap<int, int> _cooC;
    // blockBmm(bcsrA, bcscB, _cooC);
    maskedBlockBmm(bcsrA, bcsrA, bcscB, _cooC);

    t = util::toc(timer);
    std::cout << "\nBlock-BMM completed\n" << "Block-BMM time = " << t << " seconds" << std::endl;

    std::vector<std::pair<int, int>> _vecC;


    timer = util::tic();

    for (const auto& x : _cooC) {
      _vecC.push_back(std::pair<int, int> (x.first, x.second));
    }
    std::sort(_vecC.begin(), _vecC.end());

    t = util::toc(timer);
    std::cout << "\nVector processing time = " << t << " seconds" << std::endl;

    // prt::vec(_vecC);

    /* ----------------------------- gather results ----------------------------- */

    /* -------------------------------------------------------------------------- */
    /*                                    TODO                                    */
    /* -------------------------------------------------------------------------- */

    // gather _vecC from each process to process 0

    /* ------------------------------- construct C ------------------------------ */

    /* -------------------------------------------------------------------------- */
    /*                                    TODO                                    */
    /* -------------------------------------------------------------------------- */

    // having the local result of each process, add row offsets to each result and 
    // add them to the result C (C can be a vector of pairs to be sorted easier and used by the tester)

    /* ------------------------------ MPI finalize ------------------------------ */

    MPI_Finalize();

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);

    // /* ------------------------------ check result ------------------------------ */

    // if (util::checkRes(graph, vecC)) {
    //   std::cout << "\nTest passed\n";
    // }
    // else {
    //   std::cout << "\nTest failed\n";
    // }

  return 0;
}