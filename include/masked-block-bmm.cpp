/* -------------------------------------------------------------------------- */
/*                            masked-block-bmm.cpp                            */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <bits/stdc++.h>

/* ---------------------------- masked block-bmm ---------------------------- */

void maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B, std::multimap<int, int> &C)
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

    int blocksPerRow = A.n / A.b;

    // high level matrix multiplication
    for (int blockRowF = 0; blockRowF < blocksPerRow; blockRowF++) {
        for (int indF = F.HL_bRowPtr[blockRowF]; indF < F.HL_bRowPtr[blockRowF + 1]; indF++) {

            int blockColF = F.HL_bColInd[indF];
            std::multimap <int, int> _C; // block of matrix C

            maskedBlockRowColMult(blockRowF, blockColF, F, A, B, _C);
            util::addCooBlockToMatrix(C, blockRowF, blockColF, A.b, _C);
        }
    }
}

void maskedBlockRowColMult(int blockRowF, int blockColF, bcsr &F, bcsr &A, bcsc &B, std::multimap<int, int> &_C)
{   
    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRow = A.n / A.b;
    int bIndF = blockRowF * blocksPerRow + blockColF;

    int LL_rowPtrOffsetF, LL_colIndOffsetF;
    int LL_rowPtrOffsetA, LL_colIndOffsetA;
    int LL_colPtrOffsetB, LL_rowIndOffsetB;

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

            maskedBbm(F, A, B, LL_rowPtrOffsetF, LL_colIndOffsetF, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB, _C);
            
            ptr1++;
            ptr2++;
        }
    }
}

void maskedBbm( bcsr &F,
                bcsr &A,
                bcsc &B,
                int LL_rowPtrOffsetF,
                int LL_colIndOffsetF,
                int LL_rowPtrOffsetA,
                int LL_colIndOffsetA,
                int LL_colPtrOffsetB,
                int LL_rowIndOffsetB,
                std::multimap <int, int> &_C )
// masked boolean block-block multiplication
{   
    struct timeval t0 = util::tic();
    bool isExistingElement = false;

    for (int _rowF = 0; _rowF < F.b; _rowF++) {
        for (int _indF = F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF]; _indF < F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF + 1]; _indF++) {
            
            int _colF = F.LL_bColInd[_indF + LL_colIndOffsetF];

            // check if index is already true
            // auto it = _C.find(_rowF);
            // if (it != _C.end()) {
            //     if (it->second == _colF) {
            //         continue;
            //     }
            // }

            if (rowColMult(_rowF, _colF, A, B, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB )) {

                auto iter = _C.find(_rowF);
                if (iter != _C.end()) {

                    auto itr1 = _C.lower_bound(_rowF);
                    auto itr2 = _C.upper_bound(_rowF);
                    
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
                _C.insert(std::pair <int, int> (_rowF, _colF));
            }
        }
    }
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