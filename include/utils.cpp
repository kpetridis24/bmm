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
        for (int i = 0; i < len; i++)
            std::cout << arr[i] << " ";
        std::cout << std::endl;
    }


    void mat(int **mat, int rows, int cols)
    {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << mat[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void vec(std::vector<std::pair<int, int>> vec)
    {
        for (auto x : vec)
        std::cout << x.first << " " << x.second << std::endl;
    }

    void csrMat(csr &M)
    {
        std::cout << "\nrowPtr:";
        prt::arr(M.rowPtr, M.m + 1);
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

    void addCooBlockToMatrix(std::multimap<int, int> &mapM, int blockRow, int blockCol, int b, std::multimap<int, int> &_mapM)
    {
        int rowOffset = blockRow * b;
        int colOffset = blockCol * b;

        for (const auto& x : _mapM) {
            mapM.insert(std::pair <int, int> (x.first + rowOffset, x.second + colOffset));
        }
    }

    void removeCooRowOffsets(coo &M, int offset)
    {
        for (int i = 0; i < M.nnz; i++) {
            M.row[i] -= offset;
        }
    }

    void addCooRowOffsets(std::vector<std::pair <int, int>> &vecCooM, int *rowsM, int *colsM, int offset)
    {
        for (int i = 0; i < vecCooM.size(); i++) {
            rowsM[i] = vecCooM[i].first + offset;
            colsM[i] = vecCooM[i].second;
        }
    }

    void initCsr(csr &M, int m, int n, int nnz)
    // initialize CSR matrix
    {
        M.rowPtr = new int[m + 1]();
        M.colInd = new int[nnz]();
        M.m = m;
        M.n = n;
        M.nnz = nnz;
    }

    void initCsc(csc &M, int m, int n, int nnz)
    // initialize CSR matrix
    {
        M.colPtr = new int[n + 1]();
        M.rowInd = new int[nnz]();
        M.m = m;
        M.n = n;
        M.nnz = nnz;
    }

    void initCoo(coo &M, int m, int n, int nnz)
    // initialize COO matrix
    {
        M.row = new int[nnz]();
        M.col = new int[nnz]();
        M.m = m;
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

    bool checkRes(int graphInd, std::vector <std::pair <int, int>> &vecC)
    {
        std::string checkGraph;

        switch(graphInd) {
            case 0:
                checkGraph = "s6.mtx";
                break;
            case 1:
                checkGraph = "s12.mtx";
                break;
            case 2:
                checkGraph = "com-Youtube.mtx";
                break;
            case 3:
                checkGraph = "belgium_osm.mtx";
                break;
            case 4:
                checkGraph = "dblp-2010.mtx";

                break;
            case 5:
                checkGraph = "as-Skitter.mtx";
                break;
            default:
                exit(1);
        }

        // read correct result
        int checkN;
        int checkNnz;

        std::string checkFile = "graphs/bmm-res/bmm_res_" + checkGraph;

        readMtxValues(checkFile, checkN, checkNnz);
        coo checkM;
        util::initCoo(checkM, checkN, checkN, checkNnz);

        openMtxFile(checkFile, checkM.col, checkM.row, checkM.n, checkM.nnz);

        // compare results
        bool pass = true;
        if (checkM.nnz != vecC.size()) {
            return false;
        }
        else {
            for (int i = 0; i < vecC.size(); i++) {
                if (checkM.row[i] != vecC[i].first || checkM.col[i] != vecC[i].second) {
                    return false;
                }
            }
        }
        util::delCoo(checkM);
        return true;
    }
}