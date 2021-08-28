/* -------------------------------------------------------------------------- */
/*                                  utils.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>
#include <bits/stdc++.h>



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
        std::cout << "\ncolPtr:";
        prt::arr(M.colPtr, M.n + 1);
        std::cout << "rowInd:";
        prt::arr(M.rowInd, M.nnz);
    }

    void cooMat(coo &M)
    {
        std::cout << "\nrow:";
        prt::arr(M.row, M.nnz);
        std::cout << "col:";
        prt::arr(M.col, M.nnz);
    }

    void map(std::multimap <int, int> m)
    {
        for (const auto& x : m) {
            std::cout << x.first << ": " << x.second << "\n";
        }
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

    void insertCooElementInIndex(int row, int col, int ind, int *M, int &sizeM)
    {
        //std::cout << "index = " << ind << std::endl;

        // // move elements right to make space for the adding element
        // std::cout<<"size="<<sizeM<<", ptr="<<ind<<std::endl;
        int i;
        for (i = sizeM; i > ind; i-=2) {
            M[i + 1] = M[i - 1];
            M[i] = M[i - 2];
        }
        // // add element and increase matrix size
        M[i] = row;
        M[i + 1] = col;
        sizeM += 2;
    }

    void addCooElement(int row, int col, int *M, int &sizeM)
    // add an element in the right position of a COO matrix M = [row1, col1, row2, col2, ...]
    {   
        //std::cout<<"("<<row<<" , "<<col<<")"<<std::endl;

        if (sizeM == 0) { // if matrix is empty, add element
            insertCooElementInIndex(row, col, 0, M, sizeM);
            return;
        }

        for(int ptr = 0; ptr <= sizeM; ptr += 2) {

            if (ptr == sizeM) {
                util::insertCooElementInIndex(row, col, ptr, M, sizeM);
            }

            if(row > M[ptr]) {
                //std::cout<<"enter1\n";
                continue;
            }
            else if(row < M[ptr]) {
                //std::cout<<"enter2\n";
                util::insertCooElementInIndex(row, col, ptr, M, sizeM);
                return;
            }
            else { // row == M[ptr] -> search index based on col
                //std::cout<<"enter3\n";
                if(col > M[ptr + 1]) {
                    continue;
                }
                else if(col < M[ptr + 1]) {
                    // Insert here
                    util::insertCooElementInIndex(row, col, ptr, M, sizeM);
                    return;
                }
                else { // element already exists
                    return;
                }   
            }   
        }
    }

    void addCooBlockToMatrix(int *M, int *_M, int blockRow, int blockCol, int b, int &sizeM, int _sizeM)
    {
        int rowOffset = blockRow * b;
        int colOffset = blockCol * b;

        for (int i = 0; i < _sizeM; i += 2) {
            M[sizeM + i] = _M[i] + rowOffset;
            M[sizeM + i + 1] = _M[i + 1] + colOffset;
        }

        sizeM += _sizeM;
    }

    // bool searchCooElement(int row, int col, int *M, int &sizeM)
    // // search an element in a COO matrix M = [row1, col1, row2, col2, ...]
    // {
    //     for (int ptr = 0; ptr < sizeM; ptr += 2) {
    //         if (row == M[ptr]) {    // row found
    //             while (row == M[ptr]) {
    //                 if (col == M[ptr + 1])  // element found, so return true
    //                     return true;
    //                 ptr += 2;
    //             }
    //             return false;   // element is not found in its the row, so return false
    //         }
    //     }
    //     return false;
    // }

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

    bool checkRes(std::string checkGraph, coo &C)
    {
        std::vector<std::pair <int, int> > vectC (C.nnz);

        for (int i = 0; i < C.nnz; i++) {
        vectC[i].first = C.row[i];
        vectC[i].second = C.col[i];
        }

        std::sort(vectC.begin(), vectC.end());

        // read correct result

        int checkN;
        int checkNnz;

        std::string checkFile = "graphs/bmm-res/bmm_res_" + checkGraph;

        readMtxValues(checkFile, checkN, checkNnz);

        // std::cout << checkN << "\t" << checkNnz << std::endl;

        coo checkM;
        util::initCoo(checkM, checkN, checkNnz);

        openMtxFile(checkFile, checkM.col, checkM.row, checkM.n, checkM.nnz);

        // prt::cooMat(checkM);

        bool pass = true;
        if (checkM.nnz != vectC.size()) {
            return false;
        }
        else {
            for (int i = 0; i < vectC.size(); i++) {
                if (checkM.row[i] != vectC[i].first || checkM.col[i] != vectC[i].second) {
                    return false;
                }
            }
        }

        // prt::cooMat(checkM);
        util::delCoo(checkM);

        return true;
    }
}