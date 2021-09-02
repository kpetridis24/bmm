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

    /* --------------------------- Initialize OpenMPI --------------------------- */

    int numProcesses, rank;
    int graphInd = 2;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    coo cooA;
    coo cooB;

    /* ----------------------- Matrix A distribution test ----------------------- */
    coo _cooA;

    distributeCooMatrix(numProcesses, rank, cooA, _cooA, graphInd);

    broadcastCooMatrix(numProcesses, rank, cooB, graphInd);

    MPI_Finalize();

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);


    /* -------------------------------------------------------------------------- */
    /*                               block BMM test                               */
    /* -------------------------------------------------------------------------- */

    // struct timeval timer;
    // double t = -1;

    // csr A;
    // csc B;
    // int n;
    // int nnz;
    // int b = 2;

    // /* ------------------------------ read matrices ----------------------------- */

    // std::string graph = read2csr(1, n, nnz, b, A, B);

    // prt::csrMat(A);

    // /* --------------------------------- block A -------------------------------- */

    // timer = util::tic();

    // bcsr blA;
    // int numBlocks = (A.m / b) * (A.m / b);
    // int LL_bRowPtrSize = numBlocks * (b + 1);
    // blA.m = A.m;
    // blA.n = A.n;
    // blA.b = b;
    // blA.LL_bRowPtr = new int[LL_bRowPtrSize]();
    // blA.LL_bColInd = new int[A.nnz]();

    // ret _ret = csr2bcsr(A, blA);

    // blA.HL_bRowPtr = _ret.ret1;
    // blA.HL_bColInd = _ret.ret2;
    // blA.nzBlockIndex = _ret.ret3;
    // blA.blockNnzCounter = _ret.ret4;

    // t = util::toc(timer);
    // std::cout << "\nBlocking A in B-CSR completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    // /* --------------------------------- block B -------------------------------- */


    // timer = util::tic();
    // int LL_bColPtrSize = numBlocks * (b + 1);

    // bcsc blB;
    // blB.m = B.m;
    // blB.n = B.n;
    // blB.b = b;

    // // init Low-Level CSC
    // blB.LL_bColPtr = new int[LL_bColPtrSize]();
    // blB.LL_bRowInd = new int[nnz]();

    // // blocking
    // _ret = csr2bcsc(B, blB);

    // blB.HL_bColPtr = _ret.ret1;
    // blB.HL_bRowInd = _ret.ret2;
    // blB.nzBlockIndex = _ret.ret3;
    // blB.blockNnzCounter = _ret.ret4;

    // t = util::toc(timer);
    // std::cout << "\nBlocking B in B-CSC completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    // /* -------------------------------- block BMM ------------------------------- */

    // timer = util::tic();

    // std::multimap<int, int> C;
    // // blockBmm(blA, blB, C);
    // maskedBlockBmm(blA, blA, blB, C);

    // t = util::toc(timer);
    // std::cout << "\nBlock-BMM completed\n" << "Block-BMM time = " << t << " seconds" << std::endl;

    // std::vector<std::pair<int, int>> vecC;

    // for (const auto& x : C) {
    //   vecC.push_back(std::pair<int, int> (x.first, x.second));
    // }
    // std::sort(vecC.begin(), vecC.end());

    // prt::vec(vecC);

    // /* ------------------------------ check result ------------------------------ */

    // if (util::checkRes(graph, vecC)) {
    //   std::cout << "\nTest passed\n";
    // }
    // else {
    //   std::cout << "\nTest failed\n";
    // }

  return 0;
}