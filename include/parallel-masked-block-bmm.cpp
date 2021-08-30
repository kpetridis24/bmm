/* -------------------------------------------------------------------------- */
/*                        parallel-masked-block-bmm.cpp                       */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

/* ---------------------------- masked block-bmm ---------------------------- */

ret2 parallelMaskedBlockBmm(bcsr &F, bcsr &A, bcsc &B)
// masked boolean matrix multiplication F.*(A*B) using blocks
{
    if (A.n != B.n || A.n != F.n) {
        std::cout << "Dimensions error\n";
        exit(1);
    }

    if (A.b != B.b || A.b != F.b) {
        std::cout << "Block size error\n";
        exit(1);
    }

    int nnzF = F.blockNnzCounter[(F.n / F.b) * (F.n / F.b)];
    int blocksPerRow = A.n / A.b;

    int *C = new int[nnzF](); 
    int sizeC = 0;

    // high level matrix multiplication
    for (int blockRowF = 0; blockRowF < blocksPerRow; blockRowF++) {    // each iteration computes a block-row of C matrix
    
        int numOfNzBlocks = F.HL_bRowPtr[blockRowF + 1] - F.HL_bRowPtr[blockRowF];

        ret2 **bRowC = new ret2*[numOfNzBlocks];
        int startInd = F.HL_bRowPtr[blockRowF];

        cilk_for (int indF = F.HL_bRowPtr[blockRowF]; indF < F.HL_bRowPtr[blockRowF + 1]; indF++) {     // each iteration computes a specific block of the current block-row, will be parallelized
            // if we find how many iterations are executed (numOfNzrBlocks) we can make an array to store the product blocks _C and when all threads finish use addCooBlockToMatrix to add them to C
            int blockColF = F.HL_bColInd[indF];
            std::multimap <int, int> map;

            bRowC[indF - startInd] = maskedBlockRowColMult(blockRowF, blockColF, F, A, B, map);
        }

        for (int i = 0; i < numOfNzBlocks; i++) {
            int blockColF = F.HL_bColInd[F.HL_bRowPtr[blockRowF] + i];
            util::addCooBlockToMatrix(C, bRowC[i]->M, blockRowF, blockColF, A.b, sizeC, bRowC[i]->sizeM);
            delete[] bRowC[i]->M;
            delete bRowC[i];
        }

        delete[] bRowC;
    }

    ret2 ret;
    ret.M = C;
    ret.sizeM = sizeC;

    return ret;
}