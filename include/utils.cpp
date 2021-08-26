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
            std::cout << arr[i] << " ";
        std::cout << std::endl;
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

    struct timeval tic()
    {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
    }
    
    static double toc(struct timeval begin)
    {
    struct timeval end;
    gettimeofday(&end, NULL);
    double stime = ((double) (end.tv_sec - begin.tv_sec) * 1000 ) +
        ((double) (end.tv_usec - begin.tv_usec) / 1000 );
    stime = stime / 1000;
    return(stime);
    }

    void blockOffsets(int blockInd, int *nzBlockIndex, int *blockNnzCounter, int b, int &LL_rowPtrOffset, int &LL_colIndOffset)
    /* -------------------------------------------------------------------------- */
    /*             find the offsets of a specific block in the LL-CSR             */
    /* -------------------------------------------------------------------------- */
    {
        int newBlockInd =  nzBlockIndex[blockInd];
        LL_rowPtrOffset = newBlockInd * (b + 1);
        LL_colIndOffset = blockNnzCounter[blockInd];

    /* -------------------------------------------------------------------------- */
    /*                                    TODO                                    */
    /* -------------------------------------------------------------------------- */
    /* -------------------------------------------------------------------------- */
    /*          pop the nz blocks of blockNnzCounter and use newBlockInd          */
    /* -------------------------------------------------------------------------- */
    }

    void addCooElement(int row, int col, int *M, int &sizeM)
    // add an element in the right position of a COO matrix M = [row1, col1, row2, col2, ...]
    {
        if (sizeM == 0) { // if matrix is empty, add element
            M[0] = row;
            M[1] = col;
            sizeM += 2;
            return;
        }

        int ptr = 0;

        while(row > M[ptr])
            ptr += 2;

        if(row < M[ptr]) { // add element here
            //std::cout<<"enter1\n";
            int i;
            // move elements right to make space for the adding element
            for(i = sizeM; i > ptr; i-=2) {
                M[i + 1] = M[i - 1];
                M[i] = M[i - 2];
            }
            // add element and increase matrix size
            M[i] = row;
            M[i + 1] = col;
            sizeM += 2;
        }
        else { // row == M[ptr] so check the cols
            //std::cout<<"enter2\n";
            while(row == M[ptr]) {
                if (col > M[ptr + 1]) {
                    //std::cout<<"enter21\n";
                    ptr += 2;

                    if (ptr == sizeM) { // end of matrix reached so add element here
                        M[ptr] = row;
                        M[ptr + 1] = col;
                        sizeM += 2;
                        break;
                    }

                    if (row < M[ptr]) { // row changed so add element here
                        int i;
                        // move elements right to make space for the adding element
                        for(i = sizeM; i > ptr; i-=2) {
                            M[i + 1] = M[i - 1];
                            M[i] = M[i - 2];
                        }
                        // add element and increase matrix size
                        M[i] = row;
                        M[i + 1] = col;
                        sizeM += 2;
                        break;
                    }
                }
                else if ((col < M[ptr + 1])) { // add element here
                    //std::cout<<"enter22\n";
                    int i;
                    // move elements right to make space for the adding element
                    for(i = sizeM; i > ptr; i-=2) {
                        M[i + 1] = M[i - 1];
                        M[i] = M[i - 2];
                    }
                    // add element and increase matrix size
                    M[i] = row;
                    M[i + 1] = col;
                    sizeM += 2;
                    break;
                }
                else {
                    //std::cout<<"enter23\n";
                    break; // element already exists
                }
            }
        }   
    }

    bool searchCooElement(int row, int col, int *M, int &sizeM)
    // search an element in a COO matrix M = [row1, col1, row2, col2, ...]
    {
        for (int ptr = 0; ptr < sizeM; ptr += 2) {
            if (row == M[ptr]) {    // row found
                while (row == M[ptr]) {
                    if (col == M[ptr + 1])  // element found, so return true
                        return true;
                    ptr += 2;
                }
                return false;   // element is not found in its the row, so return false
            }
        }
        return false;
    }

    void addCooBlockToMatrix(int *M, int *_M, int blockRow, int blockCol, int b, int &sizeM, int _sizeM)
    {
        int rowOffset = blockRow * b;
        int colOffset = blockCol * b;

        for (int i = 0; i < _sizeM; i += 2) {
            util::addCooElement(_M[i] + rowOffset, _M[i + 1] + colOffset, M, sizeM);
        }
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

    void delBcsr(bcsr &M)
    // delete B-CSR matrix
    {
        delete[] M.LL_bRowPtr;
        delete[] M.LL_bColInd;
        delete[] M.HL_bRowPtr;
        delete[] M.HL_bColInd;
        delete[] M.nzBlockIndex;
        delete[] M.blockNnzCounter;
    }

    void delBcsc(bcsc &M)
    // delete B-CSC matrix
    {
        delete[] M.LL_bColPtr;
        delete[] M.LL_bRowInd;
        delete[] M.HL_bColPtr;
        delete[] M.HL_bRowInd;
        delete[] M.nzBlockIndex;
        delete[] M.blockNnzCounter;
    }
}