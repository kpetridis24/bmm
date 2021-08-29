/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <bits/stdc++.h>

/* ---------------------------- masked block-bmm ---------------------------- */

ret2 maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B)
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
    for (int blockRowF = 0; blockRowF < blocksPerRow; blockRowF++) {
        for (int indF = F.HL_bRowPtr[blockRowF]; indF < F.HL_bRowPtr[blockRowF + 1]; indF++) {

            int blockColF = F.HL_bColInd[indF];
            std::multimap <int, int> map;

            ret2 _C = maskedBlockRowColMult(blockRowF, blockColF, F, A, B, map);
            util::addCooBlockToMatrix(C, _C.M, blockRowF, blockColF, A.b, sizeC, _C.sizeM);
        }
    }

    ret2 ret;
    ret.M = C;
    ret.sizeM = sizeC;

    return ret;
}

ret2 maskedBlockRowColMult(int blockRowF, int blockColF, bcsr &F, bcsr &A, bcsc &B, std::multimap<int, int> &map)
{   
    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRow = A.n / A.b;
    int bIndF = blockRowF * blocksPerRow + blockColF;
    int _nnzF = F.blockNnzCounter[bIndF + 1] - F.blockNnzCounter[bIndF];
    F.nnz = _nnzF;

    int LL_rowPtrOffsetF, LL_colIndOffsetF;
    int LL_rowPtrOffsetA, LL_colIndOffsetA;
    int LL_colPtrOffsetB, LL_rowIndOffsetB;

    int *_C = new int[2 * _nnzF](); 
    int _sizeC = 0;

    int blockRowA = blockRowF;
    int blockColB = blockColF;

    while (A.HL_bRowPtr[blockRowA] + ptr1 < A.HL_bRowPtr[blockRowA + 1] && B.HL_bColPtr[blockColB] + ptr2 < B.HL_bColPtr[blockColB + 1]) {
        if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] < B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr1++;
        }
        else if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] > B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr2++;
        }
        else {
            cN = A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1]; // common neighbor index

            bIndA = blockRowA * blocksPerRow + cN;
            bIndB = blockColB * blocksPerRow + cN;

            util::blockOffsets(bIndF, F.nzBlockIndex, F.blockNnzCounter, F.b, LL_rowPtrOffsetF, LL_colIndOffsetF);
            util::blockOffsets(bIndA, A.nzBlockIndex, A.blockNnzCounter, A.b, LL_rowPtrOffsetA, LL_colIndOffsetA);
            util::blockOffsets(bIndB, B.nzBlockIndex, B.blockNnzCounter, B.b, LL_colPtrOffsetB, LL_rowIndOffsetB);

            maskedBbm(F, A, B, _C, _sizeC, LL_rowPtrOffsetF, LL_colIndOffsetF, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB, map);
            
            ptr1++;
            ptr2++;
        }
    }

    ret2 ret;
    ret.M = _C;
    ret.sizeM = _sizeC;

    return ret;
}

bool rowColMult( int rowA, int colB, 
                 bcsr &A, bcsc &B, 
                 int LL_rowPtrOffsetA,
                 int LL_colIndOffsetA,
                 int LL_colPtrOffsetB,
                 int LL_rowIndOffsetB )
// boolean inner-block row-col multiplication - multiply rowA of a block of matrix A with colB of a block of matrix B
{
    int ptr1 = 0;
    int ptr2 = 0;
    
    while (A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 < A.LL_bRowPtr[rowA + LL_rowPtrOffsetA + 1] &&
            B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 < B.LL_bColPtr[colB + LL_colPtrOffsetB + 1]) {
        if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] <
            B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
            ptr1++;
        }
        else if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] >
            B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
            ptr2++;
        }
        else {
            return true;
        }
    }

    return false;
}

void maskedBbm( bcsr &F,
                bcsr &A,
                bcsc &B,
                int *_C,
                int &_sizeC,
                int LL_rowPtrOffsetF,
                int LL_colIndOffsetF,
                int LL_rowPtrOffsetA,
                int LL_colIndOffsetA,
                int LL_colPtrOffsetB,
                int LL_rowIndOffsetB,
                std::multimap <int, int> &map )
// masked boolean block-block multiplication
{   
    struct timeval t0 = util::tic();
    bool isExistingElement = false;

    for (int _rowF = 0; _rowF < F.b; _rowF++) {
        for (int _indF = F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF]; _indF < F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF + 1]; _indF++) {
            
            int _colF = F.LL_bColInd[_indF + LL_colIndOffsetF];

            // check if index is already true
            // auto it = map.find(_rowF);
            // if (it != map.end()) {
            //     if (it->second == _colF) {
            //         continue;
            //     }
            // }

            if (rowColMult(_rowF, _colF, A, B, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB )) {

                auto iter = map.find(_rowF);
                if (iter != map.end()) {

                    auto itr1 = map.lower_bound(_rowF);
                    auto itr2 = map.upper_bound(_rowF);
                    
                    while (itr1 != itr2)
                    {   
                        if (itr1 -> first == _rowF && itr1 -> second == _colF) { // already exists.
                            isExistingElement = true;
                            break;
                        }
                        itr1++;
                    }
                 }
                if (isExistingElement) {
                    isExistingElement = false;
                    continue;
                }
                map.insert(std::pair <int, int> (_rowF, _colF));
               
                if (_sizeC <= F.nnz + 2) {
                    _C[_sizeC++] = _rowF;
                    _C[_sizeC++] = _colF;
                }
            }
        }
    }
}