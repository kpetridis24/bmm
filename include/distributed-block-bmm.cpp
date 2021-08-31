/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <headers.hpp>
#include <mpi.h>

MPI_Datatype createMpiTypeBcsc( int LL_bColPtrSize, 
                        int LL_bRowIndSize,
                        int HL_bColPtrSize,
                        int HL_bRowIndSize,
                        int nzBlockIndexSize,
                        int blockNnzCounterSize ) {

    const int nitems = 8;
    int blocklengths[8] = {LL_bColPtrSize, LL_bRowIndSize, HL_bColPtrSize, HL_bRowIndSize,
        nzBlockIndexSize, blockNnzCounterSize, 1, 1};
    MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Datatype MPI_BCSC;
    MPI_Aint offsets[8];

    offsets[0] = offsetof(bcsc, LL_bColPtr);
    offsets[1] = offsetof(bcsc, LL_bRowInd);
    offsets[2] = offsetof(bcsc, HL_bColPtr);
    offsets[3] = offsetof(bcsc, HL_bRowInd);
    offsets[4] = offsetof(bcsc, nzBlockIndex);
    offsets[5] = offsetof(bcsc, blockNnzCounter);
    offsets[6] = offsetof(bcsc, n);
    offsets[7] = offsetof(bcsc, b);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_BCSC);
    MPI_Type_commit(&MPI_BCSC);

    return MPI_BCSC;
}