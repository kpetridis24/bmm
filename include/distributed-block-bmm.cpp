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
    int *numBlocksInChunk = new int[numProcesses];
    int *offsets_LL_bRowPtr = new int[numProcesses + 1];
    int *offsets_LL_bColInd = new int[numProcesses + 1];
    int *offsets_chunkBlocks = new int[numProcesses + 1];
    int self_chunkSize_LL_bRowPtr;
    int self_chunkSize_LL_bColInd;
    int self_chunkSize_blockNnzCounter;
    int blockRow;

    if(rank == 0) {

        blockRow = 0;
        offsets_LL_bRowPtr[0] = 0;
        offsets_LL_bColInd[0] = 0;
        offsets_chunkBlocks[0] = 0;

        for (int chunkInd = 0; chunkInd < numProcesses; chunkInd++) {
            
            // find the number of blocks of the current chunk
            numBlocksInChunk[chunkInd] = 0;
            for (int _blockRow = 0; _blockRow < blockRowChunkSize; _blockRow++) {

                numBlocksInChunk[chunkInd] += M.HL_bRowPtr[blockRow + 1] - M.HL_bRowPtr[blockRow];
                offsets_chunkBlocks[chunkInd + 1] = offsets_chunkBlocks[chunkInd] + numBlocksInChunk[chunkInd];
                blockRow++;
            }

            // std::cout << numBlocksInChunk[chunkInd] << " blocks in chunk " << chunkInd << std::endl;

            chunkSizes_LL_bRowPtr[chunkInd] = numBlocksInChunk[chunkInd] * (M.b + 1);
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

        // for (int i = 0; i < (M.n / M.b) * (M.n / M.b); i ++) {
        //     if(M.blockNnzCounter[i + 1] == M.blockNnzCounter[i]) {

        //     }
        // }
    
        // prt::arr(numBlocksInChunk, numProcesses);

        // fix blockNnzCounter, add empty blocks
        int curChunk = -1;
        for (int i = 0, chunkInd = 0, cnt = 0; i < (M.n / M.b) * (M.n / M.b) + 1; i++) {

            if (M.blockNnzCounter[i] == M.blockNnzCounter[i + 1]) {
                // find chunk of index i 
                for (int j = 0; j < numProcesses; j++) {
                    if (i < offsets_chunkBlocks[j + 1]) {
                        curChunk = j;
                        break;
                    }
                }

                numBlocksInChunk[curChunk]++;
                for (int k = curChunk + 1; k < numProcesses + 1; k++) {
                    offsets_chunkBlocks[k]++;
                }
            }
        }

        // magic line
        numBlocksInChunk[numProcesses - 1]++;

        self_chunkSize_LL_bRowPtr = chunkSizes_LL_bRowPtr[0];
        // prt::arr(numBlocksInChunk, numProcesses);

        /* ----------------- send chunk sizes to the other processes ---------------- */

        for(int i = 1; i < numProcesses; i++) {

            MPI_Isend(&chunkSizes_LL_bRowPtr[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
        }

        self_chunkSize_LL_bColInd = chunkSizes_LL_bColInd[0];

        for(int i = 1; i < numProcesses; i++) {

            MPI_Isend(&chunkSizes_LL_bColInd[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &req);
        }

        self_chunkSize_blockNnzCounter = numBlocksInChunk[0];

        for(int i = 1; i < numProcesses; i++) {

            MPI_Isend(&numBlocksInChunk[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &req);
        }

    }
    else{
        /* ------------------- receive chunk sizes from process 0 ------------------- */

        MPI_Recv(&self_chunkSize_LL_bRowPtr, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        M.LL_bRowPtr = new int[self_chunkSize_LL_bRowPtr];  // Every process allocates exact ammount of memory
        // std::cout<< self_chunkSize_LL_bRowPtr <<std::endl;

        MPI_Recv(&self_chunkSize_LL_bColInd, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
        M.LL_bColInd = new int[self_chunkSize_LL_bColInd];

        MPI_Recv(&self_chunkSize_blockNnzCounter, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &stat);
        M.blockNnzCounter = new int[self_chunkSize_blockNnzCounter];
    }

    /* --------------------------- distribute matrix A -------------------------- */

    MPI_Scatterv(M.LL_bRowPtr, chunkSizes_LL_bRowPtr, offsets_LL_bRowPtr, MPI_INT, M.LL_bRowPtr,
                 self_chunkSize_LL_bRowPtr, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(M.LL_bColInd, chunkSizes_LL_bColInd, offsets_LL_bColInd, MPI_INT, M.LL_bColInd,
                 self_chunkSize_LL_bColInd, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatterv(M.blockNnzCounter, numBlocksInChunk, offsets_chunkBlocks, MPI_INT, M.blockNnzCounter,
                 self_chunkSize_blockNnzCounter, MPI_INT, 0, MPI_COMM_WORLD);

    // prt::arr(M.LL_bRowPtr, self_chunkSize_LL_bRowPtr);
    // prt::arr(M.LL_bColInd, self_chunkSize_LL_bColInd);
    // prt::arr(M.blockNnzCounter, self_chunkSize_blockNnzCounter);

    for (int i = 1; i < self_chunkSize_blockNnzCounter; i++) {
        M.blockNnzCounter[i] -= M.blockNnzCounter[0];
    }
    M.blockNnzCounter[0] = 0;

    // prt::arr(M.blockNnzCounter, self_chunkSize_blockNnzCounter);

    delete[] chunkSizes_LL_bRowPtr;
    delete[] offsets_LL_bRowPtr;
    // TODO free memory
}





// int cnt = 0, cnt2 = 0;
        // for (int i = 0; i < numProcesses; i ++) {
            
        //     int size = numBlocksInChunk[i];

        //     for (int j = 0; j < size; j++, cnt++) {
        //     //while(cnt2 < size) {

        //         if (M.blockNnzCounter[cnt + 1] == M.blockNnzCounter[cnt]) {
        //             j--;
        //             numBlocksInChunk[i] ++;
        //         }
        //         //else {
        //          //   cnt2++;
        //         //}
        //     }
        // }