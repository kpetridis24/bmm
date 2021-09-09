/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <bits/stdc++.h>

/* ---------------------------- masked block-bmm ---------------------------- */

void blockBmm(bcsr &A, bcsc &B, std::multimap<int, int> &C)
// boolean matrix multiplication (A*B) using blocks
{
    if (A.n != B.m) {
        std::cout << "Dimensions error\n";
        exit(1);
    }

    if (A.b != B.b) {
        std::cout << "Block size error\n";
        exit(1);
    }

    int numBlockRowsA = A.m / A.b;
    int blocksPerRowA = A.n / A.b;
    int numBlockRowsB = B.m / A.b;
    int blocksPerRowB = B.n / A.b;

    // high level matrix multiplication
    for (int blockRowA = 0; blockRowA < numBlockRowsA; blockRowA++) {
        for (int blockColB = 0; blockColB < blocksPerRowB; blockColB++) {

            std::multimap <int, int> _C; // block of matrix C

            blockRowColMult(blockRowA, blockColB, A, B, _C);
            util::addCooBlockToMatrix(C, blockRowA, blockColB, A.b, _C);
        }
    }
}

void blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B, std::multimap<int, int> &_C)
{   
    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRowA = A.n / A.b;

    int LL_rowPtrOffsetA, LL_colIndOffsetA;
    int LL_colPtrOffsetB, LL_rowIndOffsetB;

    while (A.HL_bRowPtr[blockRowA] + ptr1 < A.HL_bRowPtr[blockRowA + 1] && B.HL_bColPtr[blockColB] + ptr2 < B.HL_bColPtr[blockColB + 1]) {
        if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] < B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr1++;
        }
        else if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] > B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
            ptr2++;
        }
        else {
            cN = A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1]; // common neighbor index

            bIndA = blockRowA * blocksPerRowA + cN;
            bIndB = blockColB * blocksPerRowA + cN;

            util::blockOffsets(bIndA, A.nzBlockIndex, A.blockNnzCounter, A.b, LL_rowPtrOffsetA, LL_colIndOffsetA);
            util::blockOffsets(bIndB, B.nzBlockIndex, B.blockNnzCounter, B.b, LL_colPtrOffsetB, LL_rowIndOffsetB);

            bbm(A, B, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB, _C);
            
            ptr1++;
            ptr2++;
        }
    }
}

void bbm(   bcsr &A,
            bcsc &B,
            int LL_rowPtrOffsetA,
            int LL_colIndOffsetA,
            int LL_colPtrOffsetB,
            int LL_rowIndOffsetB,
            std::multimap <int, int> &_C   )
// boolean block-block multiplication
{   
    struct timeval t0 = util::tic();
    bool isExistingElement = false;

    for (int _rowA = 0; _rowA < A.b; _rowA++) {
        for (int _colB = 0; _colB < B.b; _colB++) {
            
            // check if index is already true
            // auto it = _C.find(_rowF);
            // if (it != _C.end()) {
            //     if (it->second == _colF) {
            //         continue;
            //     }
            // }

            if (rowColMult(_rowA, _colB, A, B, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB )) {

                auto iter = _C.find(_rowA);
                if (iter != _C.end()) {

                    auto itr1 = _C.lower_bound(_rowA);
                    auto itr2 = _C.upper_bound(_rowA);
                    
                    while (itr1 != itr2)
                    {   
                        if (itr1 -> first == _rowA && itr1 -> second == _colB) { // already exists.
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
                _C.insert(std::pair <int, int> (_rowA, _colB));
            }
        }
    }
}