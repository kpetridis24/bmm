/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

//#include <unordered_map>
#include <headers.hpp>
#include <bits/stdc++.h>

// void blockBmm(bcsr &A, bcsc &B)
// // boolean matrix multiplication (A*B) using blocks
// {
//     if (A.n != B.n) {
//         std::cout << "Dimensions error\n";
//         exit(1);
//     }

//     if (A.b != B.b) {
//         std::cout << "Block size error\n";
//         exit(1);
//     }

//     int blocksPerRow = A.n / A.b;

//     // high level matrix multiplication
//     for (int blockRowA = 0; blockRowA < blocksPerRow; blockRowA++) {
//         for (int blockColB = 0; blockColB < blocksPerRow; blockColB++) {
//             // if (blockRowColMult(blockRowA, blockColB, A, B)) {
//             //     // std::cout << "C(" << blockRowA << ", " << blockColB << ") is possible nzb\n";
//             // }
//             blockRowColMult(blockRowA, blockColB, A, B);
//         }
//     }
// }

// bool *blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B)
// // boolean blockRow-blockCol multiplication - multiply blockRowA of matrix A with blockColB of matrix B
// {
//     bool *_C; // = new bool[A.b * A.b]();

//     int ptr1 = 0;
//     int ptr2 = 0;
//     int bIndA;
//     int bIndB;
//     int cN;
//     int blocksPerRow = A.n / A.b;

//     while (A.HL_bRowPtr[blockRowA] + ptr1 < A.HL_bRowPtr[blockRowA + 1] && B.HL_bColPtr[blockColB] + ptr2 < B.HL_bColPtr[blockColB + 1]) {
//         if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] < B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
//             ptr1++;
//         }
//         else if (A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1] > B.HL_bRowInd[B.HL_bColPtr[blockColB] + ptr2]) {
//             ptr2++;
//         }
//         else {
//             cN = A.HL_bColInd[A.HL_bRowPtr[blockRowA] + ptr1]; // common neighbor index

//             bIndA = blockRowA * blocksPerRow + cN;
//             bIndB = blockColB * blocksPerRow + cN;

//             int LL_rowPtrOffsetA;
//             int LL_colIndOffsetA;

//             int LL_colPtrOffsetB;
//             int LL_rowIndOffsetB;

//             util::blockOffsets(bIndA, A.nzBlockIndex, A.blockNnzCounter, A.b, LL_rowPtrOffsetA, LL_colIndOffsetA);
//             util::blockOffsets(bIndB, B.nzBlockIndex, B.blockNnzCounter, B.b, LL_colPtrOffsetB, LL_rowIndOffsetB);

//             // std::cout << "\nLow-Level Multiplication: A(" << blockRowA << ", " << cN << ")" << " * B(" << cN << ", " << blockColB << ")" << std::endl;
//             // std::cout << "LL_rowPtrOffsetA = " << LL_rowPtrOffsetA << "\t" << "LL_colPtrOffsetB = " << LL_colPtrOffsetB << std::endl;
//             // std::cout << "LL_colIndOffsetA = " << LL_colIndOffsetA << "\t" << "LL_rowIndOffsetB = "<< LL_rowIndOffsetB << std::endl;

//             bbm(A, B, _C, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB);

//             ptr1++;
//             ptr2++;
//         }
//     }

//     return _C;
// }

// void bbm(   bcsr &A,
//             bcsc &B,
//             bool *_C,
//             int LL_rowPtrOffsetA,
//             int LL_colIndOffsetA,
//             int LL_colPtrOffsetB,
//             int LL_rowIndOffsetB    )
// // boolean block-block multiplication
// {
//     csr _A;
//     _A.rowPtr = A.LL_bRowPtr + LL_rowPtrOffsetA;
//     _A.colInd = A.LL_bColInd + LL_colIndOffsetA;
//     _A.n = A.b;

//     csc _B;
//     _B.colPtr = B.LL_bColPtr + LL_colPtrOffsetB;
//     _B.rowInd = B.LL_bRowInd + LL_rowIndOffsetB;
//     _B.n = B.b;

//     for (int rowA = 0; rowA < A.b; rowA++) {
//         for (int colB = 0; colB < B.b; colB++) {
//                 if (rowColMult(rowA, colB, _A, _B)) {
//                     // _C[colB * A.b + rowA] = true;
// /* -------------------------------------------------------------------------- */
// /*                                    TODO                                    */
// /* -------------------------------------------------------------------------- */
//                 }
//         }
//     }
// }

/* -------------------------------------------------------------------------- */
/*                              masked-block-bmm                              */
/* -------------------------------------------------------------------------- */

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
    auto iter = map.begin();
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

    int *_C = new int[2 * _nnzF](); // = new bool[A.b * A.b]();
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

            // std::cout << "\nLow-Level Multiplication: A(" << blockRowA << ", " << cN << ")" << " * B(" << cN << ", " << blockColB << ")" << std::endl;
            // std::cout << "LL_rowPtrOffsetA = " << LL_rowPtrOffsetA << "\t" << "LL_colPtrOffsetB = " << LL_colPtrOffsetB << std::endl;
            // std::cout << "LL_colIndOffsetA = " << LL_colIndOffsetA << "\t" << "LL_rowIndOffsetB = "<< LL_rowIndOffsetB << std::endl;
            
            maskedBbm(F, A, B, _C, _sizeC, LL_rowPtrOffsetF, LL_colIndOffsetF, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB, map, iter);
            
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
// boolean row-col multiplication - multiply rowA of matrix A with colB of matrix B
{
    int ptr1 = 0;
    int ptr2 = 0;
    
    while (A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 < A.LL_bRowPtr[rowA + LL_rowPtrOffsetA + 1] && B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 < B.LL_bColPtr[colB + LL_colPtrOffsetB + 1]) {
        if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] < B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
            ptr1++;
        }
        else if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] > B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
            ptr2++;
        }
        else {
            return true;
        }
    }

    return false;
}

void printMap(std::multimap <int, int> m)
{
    for (const auto& x : m) {
        std::cout << x.first << ": " << x.second << "\n";
    }
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
                std::multimap <int, int> &map,
                std::multimap <int, int>::iterator it )
// masked boolean block-block multiplication
{   
    struct timeval t0 = util::tic();
    bool isExistingElement = false;

    for (int _rowF = 0; _rowF < F.b; _rowF++) {
        for (int _indF = F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF]; _indF < F.LL_bRowPtr[_rowF + LL_rowPtrOffsetF + 1]; _indF++) {
            
            int _colF = F.LL_bColInd[_indF + LL_colIndOffsetF];
                
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
                map.insert(std::pair<int, int>(_rowF, _colF));
               
                if (_sizeC <= F.nnz + 2) {
                    _C[_sizeC++] = _rowF;
                    _C[_sizeC++] = _colF;
                }
            }
        }
    }
}