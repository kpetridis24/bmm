/* -------------------------------------------------------------------------- */
/*                                   bmm.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

void bmm(csr A, csc B)
{
    if (A.n != B.n) {
        std::cout << "Dimensions error\n";
        return;
    }

    for (int rowA = 0; rowA < A.n; rowA++) {
        for (int colB = 0; colB < B.n; colB++) {
            if (rowColMult(rowA, colB, A, B))
                std::cout << rowA << "\t" << colB << std::endl;
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