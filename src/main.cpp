#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <utils.hpp>
#include <bmm.hpp>


int main(){

    int n;
    int nnz;
    int b;

    std::string graph = "s12.mtx";
    std::string file = "graphs/" + graph;

/* -------------------------------------------------------------------------- */
/*                                  read file                                 */
/* -------------------------------------------------------------------------- */

    util::readmtxvalues(file, n, nnz);

    int *coo_row = new int[nnz];
    int *coo_col = new int[nnz];

    util::openmtxfile(file, coo_col, coo_row, n, nnz);
    
/* -------------------------------------------------------------------------- */
/*                             convert COO to CSR                             */
/* -------------------------------------------------------------------------- */
    std::cout<<n<<", "<<nnz<<std::endl;
    int *rowPtr = new int[n];
    int *colInd = new int[nnz];

    util::coo2csr(rowPtr, colInd, coo_row, coo_col, nnz, n, 0);

    std::cout << "\nCSR row_ptr:\t";
    util::printArray(rowPtr, n);
    std::cout << "CSR col_ind:\t";
    util::printArray(colInd, nnz);

    delete[] coo_row;
    delete[] coo_col;
    
}