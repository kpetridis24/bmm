/* -------------------------------------------------------------------------- */
/*                                   bmm.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

// void bmm(csr &A, csc &B, coo &C)
// // boolean matrix multiplication (A*B)
// {
//     if (A.n != B.n) {
//         std::cout << "Dimensions error\n";
//         exit(1);
//     }

//     C.nnz = 0;

//     for (int rowA = 0; rowA < A.n; rowA++) {
//         for (int colB = 0; colB < B.n; colB++) {
//             if (rowColMult(rowA, colB, A, B)) {
//                 C.row[C.nnz] = rowA;
//                 C.col[C.nnz] = colB;
//                 C.nnz++;
//             }
//         }
//     }
// }

// void maskedBmm(csr &F, csr &A, csc &B, coo &C)
// // masked boolean matrix multiplication F.*(A*B)
// {
//     if (F.n != A.n || A.n != B.n) {
//         std::cout << "Dimensions error\n";
//         exit(1);
//     }
    
//     C.nnz = 0;
    
//     /* ---------------------------------- test ---------------------------------- */

//     // int sizeCounter = 0;
//     // for (int rowF = 0; rowF < F.n; rowF++) {
//     //     for (int indF = F.rowPtr[rowF]; indF < F.rowPtr[rowF + 1]; indF++) {

//     //         int colF = F.colInd[indF];
//     //         if (rowColMult(rowF, colF, A, B)) {
//     //             sizeCounter++;
//     //         }
//     //     }
//     // }

//     /* -------------------------------- end test -------------------------------- */

//     for (int rowF = 0; rowF < F.n; rowF++) {
//         for (int indF = F.rowPtr[rowF]; indF < F.rowPtr[rowF + 1]; indF++) {

//             int colF = F.colInd[indF];
//             if (rowColMult(rowF, colF, A, B)) {
//                 C.row[C.nnz] = rowF;
//                 C.col[C.nnz] = colF;
//                 C.nnz++;
//             }
//         }
//     }
// }

// bool rowColMult( int rowA, int colB, 
//                  bcsr &A, bcsc &B, 
//                  int LL_rowPtrOffsetA,
//                  int LL_colIndOffsetA,
//                  int LL_colPtrOffsetB,
//                  int LL_rowIndOffsetB )
// // boolean row-col multiplication - multiply rowA of matrix A with colB of matrix B
// {
//     int ptr1 = 0;
//     int ptr2 = 0;
    
//     while (A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 < A.LL_bRowPtr[rowA + LL_rowPtrOffsetA + 1] && B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 < B.LL_bColPtr[colB + LL_colPtrOffsetB + 1]) {
//         if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] < B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
//             ptr1++;
//         }
//         else if (A.LL_bColInd[A.LL_bRowPtr[rowA + LL_rowPtrOffsetA] + ptr1 + LL_colIndOffsetA] > B.LL_bRowInd[B.LL_bColPtr[colB + LL_colPtrOffsetB] + ptr2 + LL_rowIndOffsetB]) {
//             ptr2++;
//         }
//         else {
//             return true;
//         }
//     }

//     return false;
// }

/*
bool search(int *rowPtr, int *colInd, int col, int loc)
{
    int row = rowPtr[loc];
    int num = rowPtr[loc + 1] - rowPtr[loc];

    for(int i = 0; i < num; i++){

        if(colInd[row + i] == col)
            return true;
    }
    return false;
}

void bmm2(csr &A, csr &B, coo &C)
{
    int index, num, colA;
    C.nnz = 0;

    for(int rowA = 0; rowA < A.n; rowA++) {     // row of A = i
        for(int colB = 0; colB < A.n; colB++) {    // col of B = j

            index = A.rowPtr[rowA];
            num = A.rowPtr[rowA + 1] - A.rowPtr[rowA];

            for(int z = 0; z < num; z++){

                colA = A.colInd[index + z];
                if(search(B.rowPtr, B.colInd, colB, colA)) {
                    // std::cout << rowA <<" , "<< colB << std::endl;
                    C.row[C.nnz] = rowA;
                    C.col[C.nnz] = colB;
                    C.nnz++;
                    break; // C(rowA, colB) = 1
                }
                //else cij = 0
            }
        }
    }
}

void maskedBmm2(csr &F, csr &A, csr &B, coo &C)
{
    int index, num, rowA, colA, colB;
    C.nnz = 0;

    for (int rowF = 0; rowF < F.n; rowF++) {
        for (int indF = F.rowPtr[rowF]; indF < F.rowPtr[rowF + 1]; indF++) {
            rowA = rowF;
            colB = F.colInd[indF];

            index = A.rowPtr[rowA];
            num = A.rowPtr[rowA + 1] - A.rowPtr[rowA];

            for(int z = 0; z < num; z++){

                colA = A.colInd[index + z];
                if(search(B.rowPtr, B.colInd, colB, colA)) {
                    C.row[C.nnz] = rowA;
                    C.col[C.nnz] = colB;
                    C.nnz++;
                    // std::cout << rowA <<" , "<< colB << std::endl;
                    break; // C(rowA, colB) = 1
                }
                //else cij = 0
            }
        }
    }
}
*/