/* -------------------------------------------------------------------------- */
/*                                   bmm.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

void bmm(csr &A, csc &B, coo &C)
// boolean matrix multiplication (A*B)
{
    if (A.n != B.n) {
        std::cout << "Dimensions error\n";
        return;
    }

    C.nnz = 0;

    for (int rowA = 0; rowA < A.n; rowA++) {
        for (int colB = 0; colB < B.n; colB++) {
            if (rowColMult(rowA, colB, A, B)) {
                std::cout << rowA << "\t" << colB << std::endl;
                C.row[C.nnz] = rowA;
                C.col[C.nnz] = colB;
                C.nnz++;
            }
        }
    }
}

void maskedBmm(csr &F, csr &A, csc &B, coo &C)
// masked boolean matrix multiplication F.*(A*B)
{
    if (F.n != A.n || A.n != B.n) {
        std::cout << "Dimensions error\n";
        return;
    }

    C.nnz = 0;

    for (int rowF = 0; rowF < F.n; rowF++) {
        for (int indF = F.rowPtr[rowF]; indF < F.rowPtr[rowF + 1]; indF++) {
            int colF = F.colInd[indF];
            if (rowColMult(rowF, colF, A, B)) {
                std::cout << rowF << "\t" << colF << std::endl;
                C.row[C.nnz] = rowF;
                C.col[C.nnz] = colF;
                C.nnz++;
            }
        }
    }
}

bool rowColMult(int rowA, int colB, csr A, csc B)
// boolean row-col multiplication - multiply rowA of matrix A with colB of matrix B
{
    int ptr1 = 0;
    int ptr2 = 0;

    while (A.rowPtr[rowA] + ptr1 < A.rowPtr[rowA + 1] && B.colPtr[colB] + ptr2 < B.colPtr[colB + 1]) {
        if (A.colInd[A.rowPtr[rowA] + ptr1] < B.rowInd[B.colPtr[colB] + ptr2]) {
            ptr1++;
        }
        else if (A.colInd[A.rowPtr[rowA] + ptr1] > B.rowInd[B.colPtr[colB] + ptr2]) {
            ptr2++;
        }
        else {
            return true;
        }
    }

    return false;
}