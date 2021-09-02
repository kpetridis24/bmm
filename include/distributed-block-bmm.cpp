/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <headers.hpp>
#include <mpi.h>

void distributeCooMatrix(int numProcesses, int rank, coo &M, coo &_M, int graphInd)
{
    MPI_Request req;
    MPI_Status stat;

    int n;
    int nnz;
    int b;
    int numBlockRows;
    int bRowsPerChunk;
    int rowsPerChunk; // = _m
    int currentElement;
    int selfChunkSize;
    int chunkInd;
    struct timeval timer;
    double t = -1;

    int *chunkSizes = new int[numProcesses]();
    int *chunkOffsets = new int[numProcesses]();

    /* -------------------------------- process 0 ------------------------------- */

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        read2coo(graphInd, n, nnz, b, M);
        // std::cout << "\nMatrix read successfully\nn = " << M.n << ", nnz = " << M.nnz << std::endl;
        // prt::cooMat(M);

        /* ------------------- compute num of block rows per chunk ------------------ */
        
        timer = util::tic();
        
        numBlockRows = n / b;

        if (numBlockRows % numProcesses) {
            std::cout << "numBlockRows % numProcesses != 0 \n";
            exit(1);
        }

        bRowsPerChunk = numBlockRows / numProcesses;
        rowsPerChunk = bRowsPerChunk * b;

        // std::cout << "b = " << b << std::endl;
        // std::cout << "numProcesses = " << numProcesses << std::endl;
        // std::cout << "numBlockRows = " << numBlockRows << std::endl;
        // std::cout << "bRowsPerChunk = " << bRowsPerChunk << std::endl;

        /* -------------------- compute num of elements per chunk ------------------- */

        for (int cnt = 0, elementCount = 0; cnt < M.nnz; cnt++, elementCount++) {

            currentElement = M.row[cnt];
            chunkInd = currentElement / rowsPerChunk;
            // std::cout << "chunkId=" << chunkInd << std::endl;

            chunkSizes[chunkInd]++;
            chunkOffsets[chunkInd + 1] = chunkOffsets[chunkInd] + chunkSizes[chunkInd];
        }   

        selfChunkSize = chunkSizes[0];

        // prt::arr(chunkSizes, numProcesses);
        // prt::arr(chunkDisplacements, numProcesses);

        /* ---------- send array sizes to allocate memory on all processes ---------- */
        
        for (int p = 1; p < numProcesses; p++) 
            MPI_Send(&chunkSizes[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD);
        
    }
    else {
        /* --------------------------- receive array sizes -------------------------- */

        MPI_Recv(&selfChunkSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

    }

    /* --------------------------- broadcast variables -------------------------- */

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowsPerChunk, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // std::cout << "chunk " << rank << "\t_m = " << rowsPerChunk << "\tn = " << n << "\tb = " << b << std::endl;

    /* ---------------------------- distribute matrix --------------------------- */

    util::initCoo(_M, rowsPerChunk, n, selfChunkSize);

    MPI_Scatterv(M.row, chunkSizes, chunkOffsets, MPI_INT, _M.row,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Scatterv(M.col, chunkSizes, chunkOffsets, MPI_INT, _M.col,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);

    t = util::toc(timer);
    std::cout << "Distribution of matrix time = " << t << std::endl;

    if (rank == 0) {
        util::delCoo(M);
    }

    // prt::cooMat(_M);
}

void broadcastCooMatrix(int numProcesses, int rank, coo &M, int graphInd)
{
    int n;
    int nnz;
    int b;
    struct timeval timer;
    double t = -1;

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        read2coo(graphInd, n, nnz, b, M);
        std::cout << "\nMatrix read successfully\nn = " << M.n << ", nnz = " << M.nnz << std::endl;
        // prt::cooMat(M);
    }

    timer = util::tic();

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank != 0) 
        util::initCoo(M, n, n, nnz);

    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.row, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.col, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);

    t = util::toc(timer);
    std::cout << "Broadcast of matrix time = " << t << std::endl;

    // prt::cooMat(M);
}