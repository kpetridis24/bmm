/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <headers.hpp>
#include <mpi.h>

void distributeBcsrMatrix(int numProcesses, int numBlockRows, int rank, bcsr &M) {

    MPI_Status stat;
    MPI_Request req;

    int blockRowChunkSize = numBlockRows / numProcesses;
    int *chunkSizes_LL_bRowPtr = new int[numProcesses];
    int *chunkSizes_LL_bColInd = new int[numProcesses];
    int *offsets_LL_bRowPtr = new int[numProcesses + 1];
    int *offsets_LL_bColInd = new int[numProcesses + 1];
    int self_chunkSize_LL_bRowPtr;
    int self_chunkSize_LL_bColInd;
    int blockRow;

    if(rank == 0) {

        blockRow = 0;
        offsets_LL_bRowPtr[0] = 0;
        offsets_LL_bColInd[0] = 0;

        for (int chunkInd = 0; chunkInd < numProcesses; chunkInd++) {
            
            // find the number of blocks of the current chunk
            int numBlocksInChunk = 0;
            for (int _blockRow = 0; _blockRow < blockRowChunkSize; _blockRow++) {

                numBlocksInChunk += M.HL_bRowPtr[blockRow + 1] - M.HL_bRowPtr[blockRow];
                blockRow++;
            }

            // std::cout << numBlocksInChunk << " blocks in chunk " << chunkInd << std::endl;

            chunkSizes_LL_bRowPtr[chunkInd] = numBlocksInChunk * (M.b + 1);
            offsets_LL_bRowPtr[chunkInd + 1] = offsets_LL_bRowPtr[chunkInd] + chunkSizes_LL_bRowPtr[chunkInd];
            
            // std::cout << "chunk " << chunkInd << " offsets_LL_bRowPtr = " << offsets_LL_bRowPtr[chunkInd] << std::endl;
            // std::cout << "chunk " << chunkInd << " size = " << chunkSizes_LL_bRowPtr[chunkInd] << std::endl;

            // find the number of elements (of LL_bColInd) of the current chunk
            int numElemInChunk = 0;
            for (int ind = offsets_LL_bRowPtr[chunkInd] + M.b; ind < offsets_LL_bRowPtr[chunkInd + 1]; ind+=(M.b + 1)) {

                numElemInChunk += M.LL_bRowPtr[ind];
            }

            chunkSizes_LL_bColInd[chunkInd] = numElemInChunk;
            offsets_LL_bColInd[chunkInd + 1] = offsets_LL_bColInd[chunkInd] + chunkSizes_LL_bColInd[chunkInd];
            
            // std::cout << "chunk " << chunkInd << " offsets_LL_bColInd = " << offsets_LL_bColInd[chunkInd] << std::endl;
            // std::cout << "chunk " << chunkInd << " elements = " << chunkSizes_LL_bColInd[chunkInd] << std::endl;
        }

        self_chunkSize_LL_bRowPtr = chunkSizes_LL_bRowPtr[0];

        for(int i = 1; i < numProcesses; i++) {

            MPI_Isend(&chunkSizes_LL_bRowPtr[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
        }

        self_chunkSize_LL_bColInd = chunkSizes_LL_bColInd[0];

        for(int i = 1; i < numProcesses; i++) {

            MPI_Isend(&chunkSizes_LL_bColInd[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
        }

    }
    else{
        MPI_Recv(&self_chunkSize_LL_bRowPtr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        M.LL_bRowPtr = new int[self_chunkSize_LL_bRowPtr];  // Every process allocates exact ammount of memory
        // std::cout<< self_chunkSize_LL_bRowPtr <<std::endl;

        MPI_Recv(&self_chunkSize_LL_bColInd, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        M.LL_bColInd = new int[self_chunkSize_LL_bColInd];
    }

    MPI_Scatterv(M.LL_bRowPtr, chunkSizes_LL_bRowPtr, offsets_LL_bRowPtr, MPI_INT, M.LL_bRowPtr,
                 self_chunkSize_LL_bRowPtr, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(M.LL_bColInd, chunkSizes_LL_bColInd, offsets_LL_bColInd, MPI_INT, M.LL_bColInd,
                 self_chunkSize_LL_bColInd, MPI_INT, 0, MPI_COMM_WORLD);

    // prt::arr(M.LL_bRowPtr, self_chunkSize_LL_bRowPtr);
    // prt::arr(M.LL_bColInd, self_chunkSize_LL_bColInd);

    delete[] chunkSizes_LL_bRowPtr;
    delete[] offsets_LL_bRowPtr;
}
