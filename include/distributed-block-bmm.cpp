/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <headers.hpp>
#include <mpi.h>

void distributeBcsrMatrix(int numProcesses, int numBlockRows, int rank, bcsr &M) {

    MPI_Status stat;

    int blockRowChunkSize = numBlockRows / numProcesses, prevBlock, nextBlock, nnz = 0, LL_BlockOffset,
        new_LL_bRowPtrSize = 0, self_LL_bRowPtrSize;
    
    int LL_bRowPtrSizes[numProcesses], k = 0;
    int LL_bRowPtrDisplacements[numProcesses], k2 = 0;

    if(rank == 0) {

        LL_bRowPtrDisplacements[k2++] = 0;

        for(int i = 0, i2 = 1; i < numBlockRows; i++, i2++) {

            prevBlock = M.HL_bRowPtr[i];
            nextBlock = M.HL_bRowPtr[i + 1];
            nnz += nextBlock - prevBlock;
            // std::cout<<"i2="<<i2<<std::endl;

            if(i2 >= blockRowChunkSize) {
                
                LL_BlockOffset = M.HL_bRowPtr[i + 1] * (M.b + 1);
                new_LL_bRowPtrSize = M.HL_bRowPtr[i + 1] * (M.b + 1) - new_LL_bRowPtrSize;
                LL_bRowPtrSizes[k++] = new_LL_bRowPtrSize;
                LL_bRowPtrDisplacements[k2++] = LL_BlockOffset;
                // std::cout<<LL_BlockOffset<<std::endl;
                // std::cout<<"new="<<new_LL_bRowPtrSize<<std::endl;
                nnz = 0;
                i2 = 0;
                /* -------------------------------------------------------------------------- */
                /*    TODO: cut LL/HL_bRowPtr until now and send it to process, next values   */
                /*      of HL_bRowPtr must be updated, because cumulative sum has changed     */
                /* -------------------------------------------------------------------------- */

            }
        }
        
        self_LL_bRowPtrSize = LL_bRowPtrSizes[0];

        for(int i = 1; i < numProcesses; i++) {

            MPI_Send(&LL_bRowPtrSizes[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    else{
        
        MPI_Recv(&self_LL_bRowPtrSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        M.LL_bRowPtr = new int[self_LL_bRowPtrSize];  //Every process allocates exact ammount of memory
        // std::cout<<self_LL_bRowPtrSize<<std::endl;

    }
        
    MPI_Scatterv(M.LL_bRowPtr, LL_bRowPtrSizes, LL_bRowPtrDisplacements, MPI_INT, M.LL_bRowPtr,
                 self_LL_bRowPtrSize, MPI_INT, 0, MPI_COMM_WORLD);

    // prt::arr(M.LL_bRowPtr, self_LL_bRowPtrSize);
}
