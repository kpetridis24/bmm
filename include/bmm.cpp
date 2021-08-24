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
                C.row[C.nnz] = rowA;
                C.col[C.nnz] = colB;
                C.nnz++;
            }
        }
    }
}

/* -------------------------------------------------------------------------- */
/*                                    test                                    */
/* -------------------------------------------------------------------------- */
/*
bool commonNeighbors2(int rowA, int colB, csr &A, csr &B)
// commonNeighbors function
// Compare two rows (row1, row2) and return the number of common neighbors with O(k + l), if row1 has k neighbors and row2 l.
{
    int nb = 0;
    int ptr1 = 0;
    int ptr2 = 0;

    while(row_ptr[row1] + ptr1 < row_ptr[row1 + 1] && row_ptr[row2] + ptr2 < row_ptr[row2 + 1]) {
        if(col_ind[row_ptr[row1] + ptr1] < col_ind[row_ptr[row2] + ptr2]) {
            ptr1++;
        }
        else if(col_ind[row_ptr[row1] + ptr1] > col_ind[row_ptr[row2] + ptr2]) {
            ptr2++;
        }
        else {
            return true;
            ptr1++;
            ptr2++;
        }
    }

    return false;
}
*/
bool search(int *rowPtr, int *colInd, int col, int loc){

    bool check = false;
    int row = rowPtr[loc];
    int num = rowPtr[loc + 1] - rowPtr[loc];

    for(int i = 0; i < num; i++){
        if(colInd[row + i] == col){

            check = true;
            break;
        }
    }
    return check;
}

void bmm2(csr &A, csr &B)
{
    int row, num, loc;

    for(int rowA = 0; rowA < A.n; rowA++){
        for(int ind = 0; ind < A.n; ind++){

            row = A.rowPtr[rowA];
            num = A.rowPtr[rowA + 1] - A.rowPtr[rowA];

            for(int z = 0; z < num; z++){

                loc = A.colInd[row + z];
                if(search(B.rowPtr, B.colInd, ind, loc)){
                    std::cout << rowA <<" , "<< ind << std::endl;
                    break; //cij = 1
                }
                //else cij = 0
            }
        }
    }
}

/* -------------------------------------------------------------------------- */
/*                                  end test                                  */
/* -------------------------------------------------------------------------- */

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