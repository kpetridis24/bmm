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
    // find the offsets of a specific block in the LL-CSR
    {
        int newBlockInd =  nzBlockIndex[blockInd];
        LL_rowPtrOffset = newBlockInd * (b + 1);
        LL_colIndOffset = blockNnzCounter[blockInd];
    }

    void addCooBlockToMatrix(int *M, int blockRow, int blockCol, int b, int &sizeM, std::multimap<int, int> &_mapC)
    {
        int rowOffset = blockRow * b;
        int colOffset = blockCol * b;

        // for (int i = 0; i < _sizeM; i += 2) {
        //     M[sizeM + i] = _M[i] + rowOffset;
        //     M[sizeM + i + 1] = _M[i + 1] + colOffset;
        // }

        int i = 0;
        for (const auto& x : _mapC) {            
            M[sizeM + i] = x.first + rowOffset;
            M[sizeM + i + 1] = x.second + colOffset;
            i += 2;
        }

        sizeM += i;
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
    // delete CSC matrix
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
        coo checkM;
        util::initCoo(checkM, checkN, checkNnz);

        openMtxFile(checkFile, checkM.col, checkM.row, checkM.n, checkM.nnz);

        // compare results
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
        util::delCoo(checkM);
        return true;
    }
}