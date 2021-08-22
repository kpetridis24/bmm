/* -------------------------------------------------------------------------- */
/*                                  utils.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

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

    void csrMat(csr &M)
    {
        std::cout << "\nrowPtr:";
        prt::arr(M.rowPtr, M.n + 1);
        std::cout << "colInd:";
        prt::arr(M.colInd, M.nnz);
    }

    void cscMat(csc &M)
    {
        std::cout << "\nrowPtr:";
        prt::arr(M.colPtr, M.n + 1);
        std::cout << "colInd:";
        prt::arr(M.rowInd, M.nnz);
    }

    void cooMat(coo &M)
    {
        std::cout << "\nrow:";
        prt::arr(M.row, M.nnz);
        std::cout << "col:";
        prt::arr(M.col, M.nnz);
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

    void initCsr(csr &M, int n, int nnz)
    // initialize CSR matrix
    {
        M.rowPtr = new int[n + 1]();
        M.colInd = new int[nnz]();
        M.n = n;
        M.nnz = nnz;
    }

    void initCsc(csc &M, int n, int nnz)
    // initialize CSR matrix
    {
        M.colPtr = new int[n + 1]();
        M.rowInd = new int[nnz]();
        M.n = n;
        M.nnz = nnz;
    }

    void initCoo(coo &M, int n, int nnz)
    // initialize COO matrix
    {
        M.row = new int[nnz]();
        M.col = new int[nnz]();
        M.n = n;
        M.nnz = nnz;
    }

    void delCsr(csr &M)
    // delete CSR matrix
    {
        delete[] M.rowPtr;
        delete[] M.colInd;
    }

    void delCsc(csc &M)
    // delete CSR matrix
    {
        delete[] M.colPtr;
        delete[] M.rowInd;
    }

    void delCoo(coo &M)
    // delete COO matrix
    {
        delete[] M.row;
        delete[] M.col;
    }
}