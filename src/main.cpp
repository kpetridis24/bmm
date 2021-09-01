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
    /*                        OpenMPI bmm distribution test                       */
    /* -------------------------------------------------------------------------- */

    /* --------------------------- Initialize OpenMPI --------------------------- */

    int numProcesses, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request req;
    MPI_Status stat;

    // BCSR array sizes
    int LL_bRowPtrSize, LL_bColIndSize, HL_bRowPtrSize, HL_bColIndSize;
    // BCSC array sizes
    int LL_bColPtrSize, LL_bRowIndSize, HL_bColPtrSize, HL_bRowIndSize, nzBlockIndexSize, blockNnzCounterSize;

    int rxBuffer[6];
    bcsc blB;
    bcsr blA;

    int numBlockRows, blockRowChunkSize;
    struct timeval timer;
    double t = -1;

    
    /* --------------- Only process-0 runs the blocking functions --------------- */

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        int n;
        int nnz;
        int b;
        csr A;
        csc B;

        readMtx(1, n, nnz, b, A, B);

        std::cout << "\nMatrix read successfully\nn = " << A.n << ", nnz = " << A.nnz << std::endl;

        /* --------------------------- bcsr blocking test --------------------------- */

        timer = util::tic();
        int numBlocks = (n / b) * (n / b);
        LL_bRowPtrSize = numBlocks * (b + 1);

        // bcsr blA;
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

        LL_bColIndSize = nnz;
        HL_bRowPtrSize = _ret.size1;
        HL_bColIndSize = _ret.size2;

        // prt::arr(blA.HL_bRowPtr, HL_bRowPtrSize);
        // prt::arr(blA.HL_bColInd, HL_bColIndSize);
        // prt::arr(blA.LL_bRowPtr, LL_bRowPtrSize);
        // prt::arr(blA.LL_bColInd, LL_bColIndSize);

        t = util::toc(timer);
        std::cout << "\nBlocking A in B-CSR completed\n" << "Blocking time = " << t << " seconds" << std::endl;

        /* --------------------------- bcsc blocking test --------------------------- */

        timer = util::tic();

        LL_bColPtrSize = numBlocks * (b + 1);

        // bcsc blB;
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

        LL_bRowIndSize = nnz;
        HL_bColPtrSize = _ret.size1;
        HL_bRowIndSize = _ret.size2;
        nzBlockIndexSize = _ret.size3;
        blockNnzCounterSize = _ret.size4;

        numBlockRows = _ret.size1 - 1; // num_of_block_rows
        blockRowChunkSize = numBlockRows / numProcesses;
        std::cout<<"Blockrows:\t"<<numBlockRows<<std::endl<<"Blockrows per process:\t"<<
            blockRowChunkSize<<std::endl;

        /* ---------- Send array sizes to allocate memory on all processes ---------- */

        int sendBuffer[6] = {LL_bColPtrSize, LL_bRowIndSize, HL_bColPtrSize, HL_bRowIndSize,
            nzBlockIndexSize, blockNnzCounterSize};

        for(int p = 1; p < numProcesses; p++) 
            MPI_Isend(&sendBuffer, 6, MPI_INT, p, 0, MPI_COMM_WORLD, &req);
        
        t = util::toc(timer);
        std::cout << "\nBlocking B in B-CSC completed\n" << "Blocking time = " << t << " seconds" << std::endl;
    }
    else {
        /* --------------------------- Receive array sizes -------------------------- */

        MPI_Recv(&rxBuffer, 6, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

        LL_bColPtrSize = rxBuffer[0];
        LL_bRowIndSize = rxBuffer[1];
        HL_bColPtrSize = rxBuffer[2];
        HL_bRowIndSize = rxBuffer[3];
        nzBlockIndexSize = rxBuffer[4];
        blockNnzCounterSize = rxBuffer[5];

        numBlockRows = rxBuffer[2] - 1; // num_of_block_rows
        blockRowChunkSize = numBlockRows / numProcesses;

        /* ------------------- Memory allocation on all processes ------------------- */

        blB.LL_bColPtr = new int[LL_bColPtrSize]();
        blB.LL_bRowInd = new int[LL_bRowIndSize]();
        blB.HL_bColPtr = new int[HL_bColPtrSize]();
        blB.HL_bRowInd = new int[HL_bRowIndSize]();
        blB.nzBlockIndex = new int[nzBlockIndexSize]();
        blB.blockNnzCounter = new int[blockNnzCounterSize]();
    }

    /* ------------------------ BCSC B array distribution ----------------------- */

    // MPI_Bcast(blB.LL_bColPtr, LL_bColPtrSize, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(blB.LL_bRowInd, LL_bRowIndSize, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(blB.HL_bColPtr, HL_bColPtrSize, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(blB.HL_bRowInd, HL_bRowIndSize, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(blB.nzBlockIndex, nzBlockIndexSize, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(blB.blockNnzCounter, blockNnzCounterSize, MPI_INT, 0, MPI_COMM_WORLD);
    
    // prt::arr(blB.blockNnzCounter, blockNnzCounterSize);


    /* ----------------------- Matrix A distribution test ----------------------- */

    timer = util::tic();

    distributeBcsrMatrix(numProcesses, numBlockRows, rank, blA);

    t = util::toc(timer);
    std::cout << "Distribution of matrix A time = " << t << std::endl;

    MPI_Finalize();

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);

  return 0;
}