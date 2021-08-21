// /* -------------------------------------------------------------------------- */
// /*                                   bmm.cpp                                  */
// /* -------------------------------------------------------------------------- */

// #include <iostream>
// #include <headers.hpp>

// void bmm(csr A, csr B, csr C)
// {
//     if (A.n != B.n) {
//         std::cout << "Dimensions error\n";
//         return;
//     }

//     for (int rowA = 0; rowA < A.n; rowA++) {
//         for (int indA = A.rowPtr[rowA]; indA < A.rowPtr[rowA + 1]; indA++) {
//             int rowB = A.colInd[indA];
//             commonNeighbors(rowA, rowB, A, B);
//         }
//     }
// }

// int triCount(int row_ptr[], int col_ind[], int n)
// {
//     int t = 0;
//     int *c3 = new int[n]();

//     for (int row1 = 0; row1 < n; row1++) {
//         for (int idx1 = row_ptr[row1]; idx1 < row_ptr[row1 + 1]; idx1++) {
//             //int row2 = col_ind[idx1];
//             c3[row1] += commonNeighbors(row1, col_ind[idx1], row_ptr, col_ind);
//         }
//     }

//     for (int i = 0; i < n; i++) {
//         t += c3[i];
//     }

//     delete[] c3;
    
//     return t;
// }

// bool commonNeighbors(int rowA, int colB, csr A, csc B)
// // commonNeighbors function
// // Compare two rows (row1, row2) and return true if there are common vertices.
// {
//     int nb = 0;
//     int ptr1 = 0;
//     int ptr2 = 0;

//     while(row_ptr[row1] + ptr1 < row_ptr[row1 + 1] && row_ptr[row2] + ptr2 < row_ptr[row2 + 1]) {
//         if(col_ind[row_ptr[row1] + ptr1] < col_ind[row_ptr[row2] + ptr2]) {
//             ptr1++;
//         }
//         else if(col_ind[row_ptr[row1] + ptr1] > col_ind[row_ptr[row2] + ptr2]) {
//             ptr2++;
//         }
//         else {
//             return true;
//         }
//     }

//     return false;
// }