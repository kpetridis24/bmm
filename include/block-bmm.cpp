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


void maskedBlockBmm(bcsr F, bcsr A, bcsc B)
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

bool *maskedBlockRowColMult(int blockRowF, int blockColF, bcsr &F, bcsr &A, bcsc &B)
{
    bool *_C; // = new bool[A.b * A.b]();

    int ptr1 = 0;
    int ptr2 = 0;
    int bIndA;
    int bIndB;
    int cN;
    int blocksPerRow = A.n / A.b;

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

            // TODO F indices

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

            maskedBbm(A, B, _C, LL_rowPtrOffsetA, LL_colIndOffsetA, LL_colPtrOffsetB, LL_rowIndOffsetB);

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
                }
        }
    }
}

void maskedBbm( bcsr &A,
                bcsc &B,
                bool *_C,
                int LL_rowPtrOffsetA,
                int LL_colIndOffsetA,
                int LL_colPtrOffsetB,
                int LL_rowIndOffsetB )
// masked boolean block-block multiplication
{
    // TODO
}
