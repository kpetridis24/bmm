/* -------------------------------------------------------------------------- */
/*                               triCounting.cpp                              */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

/* -------------------------------------------------------------------------- */
/*                       block tri-counting using B-CSR                       */
/* -------------------------------------------------------------------------- */

int LLcommonNeighbors(  int row1, 
                        int row2, 
                        int *row_ptr1, 
                        int *col_ind1, 
                        int *row_ptr2, 
                        int *col_ind2  )
/* -------------------------------------------------------------------------- */
/*                     common nz elements between two rows                    */
/* -------------------------------------------------------------------------- */
{
    int nb = 0;
    int ptr1 = 0;
    int ptr2 = 0;

    while(row_ptr1[row1] + ptr1 < row_ptr1[row1 + 1] && row_ptr2[row2] + ptr2 < row_ptr2[row2 + 1]) {
        if(col_ind1[row_ptr1[row1] + ptr1] < col_ind2[row_ptr2[row2] + ptr2]) {
            ptr1++;
        }
        else if(col_ind1[row_ptr1[row1] + ptr1] > col_ind2[row_ptr2[row2] + ptr2]) {
            ptr2++;
        }
        else {
            nb++;
            ptr1++;
            ptr2++;
        }
    }

    return nb;
}

int LLtriCount( int row_ptr1[], 
                int col_ind1[], 
                int row_ptr2[], 
                int col_ind2[], 
                int row_ptr3[], 
                int col_ind3[], 
                int b )
/* -------------------------------------------------------------------------- */
/*                 Triangle counting between 3 blocks A, B, C                 */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*                                   A.*B*C                                   */
/* -------------------------------------------------------------------------- */
{
    int _t = 0;

    for (int row1 = 0; row1 < b; row1++) {
        for (int i = row_ptr1[row1]; i < row_ptr1[row1 + 1]; i++) {
            int col1 = col_ind1[i];
            // common neighbors B[row1, :] & C[col1, :]
            _t += LLcommonNeighbors(row1, col1, row_ptr2, col_ind2, row_ptr3, col_ind3);
        }
    }

    return _t;
}

int HLcommonNeighbors(  int row1, 
                        int row2, 
                        int *HL_b_row_ptr, 
                        int *HL_b_col_ind, 
                        int *LL_b_row_ptr, 
                        int *LL_b_col_ind, 
                        int *nzBlockIndex, 
                        int *blockNnzCounter, 
                        int n, 
                        int b )
/* -------------------------------------------------------------------------- */
/*                   common nz blocks between two block-rows                  */
/* -------------------------------------------------------------------------- */
{
    int ptr1 = 0;
    int ptr2 = 0;
    int bInd1;
    int bInd2;
    int bInd3;
    int cN;
    int blocksPerRow = n / b;

    int _t = 0;

    while (HL_b_row_ptr[row1] + ptr1 < HL_b_row_ptr[row1 + 1] && HL_b_row_ptr[row2] + ptr2 < HL_b_row_ptr[row2 + 1]) {
        if (HL_b_col_ind[HL_b_row_ptr[row1] + ptr1] < HL_b_col_ind[HL_b_row_ptr[row2] + ptr2]) {
            ptr1++;
        }
        else if (HL_b_col_ind[HL_b_row_ptr[row1] + ptr1] > HL_b_col_ind[HL_b_row_ptr[row2] + ptr2]) {
            ptr2++;
        }
        else {

            cN = HL_b_col_ind[HL_b_row_ptr[row1] + ptr1];

            // std::cout << "Common neighbor: " << cN << std::endl;

            bInd1 = row1 * blocksPerRow + row2;
            bInd2 = row1 * blocksPerRow + cN;
            bInd3 = row2 * blocksPerRow + cN;

            int LL_row_ptr_offset_1;
            int LL_col_ind_offset_1;

            int LL_row_ptr_offset_2;
            int LL_col_ind_offset_2;

            int LL_row_ptr_offset_3;
            int LL_col_ind_offset_3;

            util::blockOffsets(bInd1, nzBlockIndex, blockNnzCounter, b, LL_row_ptr_offset_1, LL_col_ind_offset_1);
            util::blockOffsets(bInd2, nzBlockIndex, blockNnzCounter, b, LL_row_ptr_offset_2, LL_col_ind_offset_2);
            util::blockOffsets(bInd3, nzBlockIndex, blockNnzCounter, b, LL_row_ptr_offset_3, LL_col_ind_offset_3);

            // std::cout << "Low-Level Multiplication: B" << bInd1 << ".*B" << bInd2 << "*B" << bInd3 << std::endl;
            // std::cout << LL_row_ptr_offset_1 << "\t" << LL_row_ptr_offset_2 << "\t" << LL_row_ptr_offset_3 << std::endl;
            // std::cout << LL_col_ind_offset_1 << "\t" << LL_col_ind_offset_2 << "\t" << LL_col_ind_offset_3 << std::endl;

            _t += LLtriCount(   LL_b_row_ptr + LL_row_ptr_offset_1, 
                                LL_b_col_ind + LL_col_ind_offset_1, 
                                LL_b_row_ptr + LL_row_ptr_offset_2, 
                                LL_b_col_ind + LL_col_ind_offset_2, 
                                LL_b_row_ptr + LL_row_ptr_offset_3, 
                                LL_b_col_ind + LL_col_ind_offset_3, 
                                b   );

            ptr1++;
            ptr2++;
        }
    }

    return _t;
}


int bCsrTriCount(   int *LL_b_row_ptr, 
                    int *LL_b_col_ind, 
                    int *HL_b_row_ptr, 
                    int *HL_b_col_ind, 
                    int *nzBlockIndex, 
                    int *blockNnzCounter,
                    int n, 
                    int b   )
{
    int blocksPerRow = n / b;
    int t = 0;
    int *c3 = new int[blocksPerRow]();

    for (int blockRow1 = 0; blockRow1 < blocksPerRow; blockRow1++) {

        // std::cout << "\nNon-zero blocks of block-row " << blockRow1 << ":\t";
        for (int i = HL_b_row_ptr[blockRow1]; i < HL_b_row_ptr[blockRow1 + 1]; i++) {   // for every nz block of blockRow1
            // std::cout << HL_b_col_ind[i] << " ";
            
            int blockRow2 = HL_b_col_ind[i];

            c3[blockRow1] += HLcommonNeighbors( blockRow1, 
                                                blockRow2, 
                                                HL_b_row_ptr, 
                                                HL_b_col_ind, 
                                                LL_b_row_ptr, 
                                                LL_b_col_ind,
                                                nzBlockIndex, 
                                                blockNnzCounter, 
                                                n, 
                                                b );

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

/* -------------------------------------------------------------------------- */
/*                             simple tri-counting                            */
/* -------------------------------------------------------------------------- */

int commonNeighbors(int row1, int row2, int row_ptr[], int col_ind[])
// commonNeighbors function
// Compare two rows (row1, row2) and return the number of common neighbors with O(k + l), if row1 has k neighbors and row2 l.
{
    int nb = 0;
    int ptr1 = 0;
    int ptr2 = 0;

    while(row_ptr[row1] + ptr1 < row_ptr[row1 + 1] && row_ptr[row2] + ptr2 < row_ptr[row2 + 1]) {
        if(col_ind[row_ptr[row1] + ptr1] < col_ind[row_ptr[row2] + ptr2]) {
            ptr1++;
        }
        else if(col_ind[row_ptr[row1] + ptr1] > col_ind[row_ptr[row2] + ptr2]) {
            ptr2++;
        }
        else {
            nb++;
            ptr1++;
            ptr2++;
        }
    }

    return nb;
}

int triCount(int row_ptr[], int col_ind[], int n)
{
    int t = 0;
    int *c3 = new int[n]();

    for (int row1 = 0; row1 < n; row1++) {
        for (int idx1 = row_ptr[row1]; idx1 < row_ptr[row1 + 1]; idx1++) {
            //int row2 = col_ind[idx1];
            c3[row1] += commonNeighbors(row1, col_ind[idx1], row_ptr, col_ind);
        }
    }

    for (int i = 0; i < n; i++) {
        t += c3[i];
    }

    delete[] c3;
    
    return t;
}