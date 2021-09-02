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

    coo cooA;
    coo cooB;

    int n;
    int nnz;
    int b;
    int numBlockRows;
    int bRowsPerChunk;
    int rowsPerChunk;
    int _elemsInChunk;
    struct timeval timer;
    double t = -1;

    /* -------------------------------- process 0 ------------------------------- */

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        read2coo(1, n, nnz, b, cooA, cooB);

        std::cout << "\nMatrix read successfully\nn = " << cooA.n << ", nnz = " << cooA.nnz << std::endl;

        /* ------------------- compute num of block rows per chunk ------------------ */

        numBlockRows = n / b;

        if (numBlockRows % numProcesses) {
            std::cout << "numBlockRows % numProcesses != 0 \n";
            exit(1);
        }

        bRowsPerChunk = numBlockRows / numProcesses;
        rowsPerChunk = bRowsPerChunk * b;

        std::cout << "b = " << b << std::endl;
        std::cout << "numProcesses = " << numProcesses << std::endl;
        std::cout << "numBlockRows = " << numBlockRows << std::endl;
        std::cout << "bRowsPerChunk = " << bRowsPerChunk << std::endl;

        /* -------------------- compute num of elements per chunk ------------------- */

        int *elemsInChunk = new int[numProcesses]();
        int chunkStartingRow = 0;
        int curRow = cooA.row[0];
        int nextChunkOffset = 0;

        for (int i = 0; i < numProcesses; i++) {  // i iterates the processes
            // std::cout << "chunk " << i << std::endl;

            chunkStartingRow = cooA.row[nextChunkOffset];
            // std::cout << "chunkStartingRow = " << chunkStartingRow << std::endl;

            for (int j = nextChunkOffset; j < cooA.nnz;) {  // j iterates the nnz of COO matrix

                curRow = cooA.row[j];

                if (curRow - chunkStartingRow < rowsPerChunk) {
                    elemsInChunk[i]++;
                    j++;
                }
                else {
                    // std::cout << "curRow - chunkStartingRow = " << curRow - chunkStartingRow << std::endl;
                    nextChunkOffset = j;
                    break;
                }
            }
        }

        _elemsInChunk = elemsInChunk[0];

        // prt::cooMat(cooA);
        // prt::arr(elemsInChunk, numProcesses);

        /* ---------- send array sizes to allocate memory on all processes ---------- */

        for(int p = 1; p < numProcesses; p++)
            MPI_Isend(&elemsInChunk[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD, &req);

        /* ---------------------- send n and b to all processes --------------------- */
        
        for(int p = 1; p < numProcesses; p++) {
            MPI_Isend(&n, 1, MPI_INT, p, 1, MPI_COMM_WORLD, &req);
            MPI_Isend(&b, 1, MPI_INT, p, 2, MPI_COMM_WORLD, &req);
        }
    }

    /* ------------------------------ process != 0 ------------------------------ */

    else {

        /* --------------------------- receive array sizes -------------------------- */

        MPI_Recv(&_elemsInChunk, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

        // std::cout << "chunk " << rank << " _elemsInChunk = " << _elemsInChunk << std::endl;

        /* ----------------------------- receive n and b ---------------------------- */

        MPI_Recv(&n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(&b, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &stat);

        // std::cout << "chunk " << rank << " n = " << n << " b = " << b << std::endl;

        /* ------------------- Memory allocation on all processes ------------------- */

        util::initCoo(cooA, n / b, _elemsInChunk);

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

    // timer = util::tic();

    // distributeBcsrMatrix(numProcesses, numBlockRows, rank, blA);

    // t = util::toc(timer);
    // std::cout << "Distribution of matrix A time = " << t << std::endl;

    MPI_Finalize();

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);

  return 0;
}