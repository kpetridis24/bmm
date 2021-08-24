/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>

int bCsrTriCount( bcsr A )
{
    int blocksPerRow = A.n / A.b;
    int t = 0;
    int *c3 = new int[blocksPerRow]();

    for (int blockRow1 = 0; blockRow1 < blocksPerRow; blockRow1++) {

        // std::cout << "\nNon-zero blocks of block-row " << blockRow1 << ":\t";
        for (int i = A.HL_bRowPtr[blockRow1]; i < A.HL_bRowPtr[blockRow1 + 1]; i++) {   // for every nz block of blockRow1
            // std::cout << A.HL_bColInd[i] << " ";
            
            int blockRow2 = A.HL_bColInd[i];

            // c3[blockRow1] += HLcommonNeighbors( blockRow1, 
            //                                     blockRow2, 
            //                                     A.HL_bRowPtr, 
            //                                     A.HL_bColInd, 
            //                                     LL_b_row_ptr, 
            //                                     LL_b_col_ind,
            //                                     nzBlockIndex, 
            //                                     blockNnzCounter, 
            //                                     n, 
            //                                     b );

            // std::cout   << "Num of common block-neighbors between block-rows " 
            //             << blockRow1 << " and " << blockRow2 << " = " 
            //             << temp
            //             << std::endl;
        }
        // std::cout << std::endl;
    }

    for (int i = 0; i < blocksPerRow; i++) {
        t += c3[i];
    }

    delete[] c3;
    
    return t;
}

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