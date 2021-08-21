/* -------------------------------------------------------------------------- */
/*                                  utils.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <iostream>

namespace prt
{
    void arr(int *arr, int len)
    {
        std::cout << std::endl;
        for(int i = 0; i < len; i++)
            std::cout << arr[i] << "\t";
        std::cout << std::endl << std::endl;
    }


    void mat(int **mat, int rows, int cols)
    {
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                std::cout << mat[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
}

namespace util
{
    void blockOffsets(int blockInd, int *nzBlockIndex, int *blockNnzCounter, int b, int &LL_row_ptr_offset, int &LL_col_ind_offset)
    /* -------------------------------------------------------------------------- */
    /*             find the offsets of a specific block in the LL-CSR             */
    /* -------------------------------------------------------------------------- */
    {
        int newBlockInd =  nzBlockIndex[blockInd];
        LL_row_ptr_offset = newBlockInd * (b + 1);
        LL_col_ind_offset = blockNnzCounter[blockInd];

    /* -------------------------------------------------------------------------- */
    /*                                    TODO                                    */
    /* -------------------------------------------------------------------------- */
    /* -------------------------------------------------------------------------- */
    /*          pop the nz blocks of blockNnzCounter and use newBlockInd          */
    /* -------------------------------------------------------------------------- */
    }
}