/* -------------------------------------------------------------------------- */
/*                                block-bmm.cpp                               */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>

void blockBmm(bcsr &A, bcsc &B)
// boolean matrix multiplication (A*B) using blocks
{
    if (A.n != B.n) {
        std::cout << "Dimensions error\n";
        exit(1);
    }

    if (A.b != B.b) {
        std::cout << "Block size error\n";
        exit(1);
    }

    int blocksPerRow = A.n / A.b;

    // high level matrix multiplication
    for (int blockRowA = 0; blockRowA < blocksPerRow; blockRowA++) {
        for (int blockColB = 0; blockColB < blocksPerRow; blockColB++) {
            // if (blockRowColMult(blockRowA, blockColB, A, B)) {
            //     // std::cout << "C(" << blockRowA << ", " << blockColB << ") is possible nzb\n";
            // }
            blockRowColMult(blockRowA, blockColB, A, B);
        }
    }
}

bool *blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B)
// boolean blockRow-blockCol multiplication - multiply blockRowA of matrix A with blockColB of matrix B
{
    bool *_C; // = new bool[A.b * A.b]();

    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRow = A.n / A.b;

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

            int LL_rowPtrOffsetA;
            int LL_colIndOffsetA;

            int LL_colPtrOffsetB;
            int LL_rowIndOffsetB;

            util::blockOffsets(bIndA, A.nzBlockIndex, A.blockNnzCounter, A.b, LL_rowPtrOffsetA, LL_colIndOffsetA);
            util::blockOffsets(bIndB, B.nzBlockIndex, B.blockNnzCounter, B.b, LL_colPtrOffsetB, LL_rowIndOffsetB);

            // std::cout << "\nLow-Level Multiplication: A(" << blockRowA << ", " << cN << ")" << " * B(" << cN << ", " << blockColB << ")" << std::endl;
            // std::cout << "LL_rowPtrOffsetA = " << LL_rowPtrOffsetA << "\t" << "LL_colPtrOffsetB = " << LL_colPtrOffsetB << std::endl;
            // std::cout << "LL_colIndOffsetA = " << LL_colIndOffsetA << "\t" << "LL_rowIndOffsetB = "<< LL_rowIndOffsetB << std::endl;

            bbm(A, B, _C, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB);

            ptr1++;
            ptr2++;
        }
    }

    return _C;
}

void bbm(   bcsr &A,
            bcsc &B,
            bool *_C,
            int LL_rowPtrOffsetA,
            int LL_colIndOffsetA,
            int LL_colPtrOffsetB,
            int LL_rowIndOffsetB    )
// boolean block-block multiplication
{
    csr _A;
    _A.rowPtr = A.LL_bRowPtr + LL_rowPtrOffsetA;
    _A.colInd = A.LL_bColInd + LL_colIndOffsetA;
    _A.n = A.b;

    csc _B;
    _B.colPtr = B.LL_bColPtr + LL_colPtrOffsetB;
    _B.rowInd = B.LL_bRowInd + LL_rowIndOffsetB;
    _B.n = B.b;

    for (int rowA = 0; rowA < A.b; rowA++) {
        for (int colB = 0; colB < B.b; colB++) {
                if (rowColMult(rowA, colB, _A, _B)) {
                    // _C[colB * A.b + rowA] = true;
/* -------------------------------------------------------------------------- */
/*                                    TODO                                    */
/* -------------------------------------------------------------------------- */
                }
        }
    }
}

/* -------------------------------------------------------------------------- */
/*                              masked-block-bmm                              */
/* -------------------------------------------------------------------------- */

void maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B)
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
            maskedBlockRowColMult(blockRowF, blockColF, F, A, B);
        }
    }
}

int *maskedBlockRowColMult(int blockRowF, int blockColF, bcsr &F, bcsr &A, bcsc &B)
{
    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRow = A.n / A.b;
    int bIndF = blockRowF * blocksPerRow + blockColF;
    int nnzF = F.blockNnzCounter[bIndF + 1] - F.blockNnzCounter[bIndF];
    int *_C = new int[2 * nnzF]; // = new bool[A.b * A.b]();
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

            int LL_rowPtrOffsetF;
            int LL_colIndOffsetF;

            int LL_rowPtrOffsetA;
            int LL_colIndOffsetA;

            int LL_colPtrOffsetB;
            int LL_rowIndOffsetB;

            util::blockOffsets(bIndF, F.nzBlockIndex, F.blockNnzCounter, F.b, LL_rowPtrOffsetF, LL_colIndOffsetF);
            util::blockOffsets(bIndA, A.nzBlockIndex, A.blockNnzCounter, A.b, LL_rowPtrOffsetA, LL_colIndOffsetA);
            util::blockOffsets(bIndB, B.nzBlockIndex, B.blockNnzCounter, B.b, LL_colPtrOffsetB, LL_rowIndOffsetB);

            // std::cout << "\nLow-Level Multiplication: A(" << blockRowA << ", " << cN << ")" << " * B(" << cN << ", " << blockColB << ")" << std::endl;
            // std::cout << "LL_rowPtrOffsetA = " << LL_rowPtrOffsetA << "\t" << "LL_colPtrOffsetB = " << LL_colPtrOffsetB << std::endl;
            // std::cout << "LL_colIndOffsetA = " << LL_colIndOffsetA << "\t" << "LL_rowIndOffsetB = "<< LL_rowIndOffsetB << std::endl;

            maskedBbm(F, A, B, _C, _sizeC, LL_rowPtrOffsetF, LL_colIndOffsetF, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB);

            ptr1++;
            ptr2++;
        }
    }
    
    delete[] _C;
/* -------------------------------------------------------------------------- */
/*                      TODO return _C and store it in C                      */
/* -------------------------------------------------------------------------- */
    return NULL;
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
                int LL_rowIndOffsetB )
// masked boolean block-block multiplication
{
    csr _F;
    _F.rowPtr = F.LL_bRowPtr + LL_rowPtrOffsetF;
    _F.colInd = F.LL_bColInd + LL_colIndOffsetF;
    _F.n = F.b;

    csr _A;
    _A.rowPtr = A.LL_bRowPtr + LL_rowPtrOffsetA;
    _A.colInd = A.LL_bColInd + LL_colIndOffsetA;
    _A.n = A.b;

    csc _B;
    _B.colPtr = B.LL_bColPtr + LL_colPtrOffsetB;
    _B.rowInd = B.LL_bRowInd + LL_rowIndOffsetB;
    _B.n = B.b;

    for (int _rowF = 0; _rowF < F.n / F.b; _rowF++) {
        for (int _indF = _F.rowPtr[_rowF]; _indF < _F.rowPtr[_rowF + 1]; _indF++) {

            int _colF = _F.colInd[_indF];
            if (rowColMult(_rowF, _colF, _A, _B)) {
                _C[_sizeC++] = _rowF;
                _C[_sizeC++] = _colF;
            }
        }
    }
}
