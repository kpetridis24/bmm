/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>

void blockBmm(bcsr &A, bcsc &B)
// boolean matrix multiplication (A*B) using blocks
{
    if (A.n != B.n) {
        std::cout << "Dimensions error\n";
        exit(1);
    }

    if (A.b != B.b) {
        std::cout << "Block size error\n";
        exit(1);
    }

    int blocksPerRow = A.n / A.b;

    csr HL_A;
    HL_A.rowPtr = A.HL_bRowPtr;
    HL_A.colInd = A.HL_bColInd;
    csc HL_B;
    HL_B.colPtr = B.HL_bColPtr;
    HL_B.rowInd = B.HL_bRowInd;

    // high level matrix multiplication
    for (int blockRowA = 0; blockRowA < blocksPerRow; blockRowA++) {
        for (int blockColB = 0; blockColB < blocksPerRow; blockColB++) {
            if (blockRowColMult(blockRowA, blockColB, A, B)) {
                // std::cout << "C(" << blockRowA << ", " << blockColB << ") is possible nzb\n";
            }
        }
    }
}


void maskedBlockBmm(bcsr F, bcsr A, bcsc B)
// masked boolean matrix multiplication F.*(A*B) using blocks
{
/* -------------------------------------------------------------------------- */
/*                                    TODO                                    */
/* -------------------------------------------------------------------------- */
}




// boolean matrix multiplication (A*B) using blocks

// int bCsrTriCount( bcsr A )
// {
//     int blocksPerRow = A.n / A.b;
//     int t = 0;
//     int *c3 = new int[blocksPerRow]();

//     for (int blockRow1 = 0; blockRow1 < blocksPerRow; blockRow1++) {

//         // std::cout << "\nNon-zero blocks of block-row " << blockRow1 << ":\t";
//         for (int i = A.HL_bRowPtr[blockRow1]; i < A.HL_bRowPtr[blockRow1 + 1]; i++) {   // for every nz block of blockRow1
//             // std::cout << A.HL_bColInd[i] << " ";
            
//             int blockRow2 = A.HL_bColInd[i];

//             // c3[blockRow1] += HLcommonNeighbors( blockRow1, 
//             //                                     blockRow2, 
//             //                                     A.HL_bRowPtr, 
//             //                                     A.HL_bColInd, 
//             //                                     LL_b_row_ptr, 
//             //                                     LL_b_col_ind,
//             //                                     nzBlockIndex, 
//             //                                     blockNnzCounter, 
//             //                                     n, 
//             //                                     b );

//             // std::cout   << "Num of common block-neighbors between block-rows " 
//             //             << blockRow1 << " and " << blockRow2 << " = " 
//             //             << temp
//             //             << std::endl;
//         }
//         // std::cout << std::endl;
//     }

//     for (int i = 0; i < blocksPerRow; i++) {
//         t += c3[i];
//     }

//     delete[] c3;
    
//     return t;
// }

bool blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B)
// boolean blockRow-blockCol multiplication - multiply blockRowA of matrix A with blockColB of matrix B
{
    int ptr1 = 0;
    int ptr2 = 0;

    while (A.HL_bRowPtr[blockRowA] + ptr1 < A.HL_bRowPtr[blockRowA + 1] && B.HL_bColPtr[blockColB] + ptr2 < B.HL_bColPtr[blockColB + 1]) {
        if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] < B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr1++;
        }
        else if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] > B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr2++;
        }
        else {
            // TODO low-level multiplication
            return true;
        }
    }

    return false;
}