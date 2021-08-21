/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include <headers.hpp>
#include <bmm.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main()
{
    /* ------------------------------- read matrix ------------------------------ */

    int n;
    int nnz;

    std::string graph = "s6.mtx";
    std::string file = "graphs/" + graph;

    read_mtx_values(file, n, nnz);

    int *coo_row = new int[nnz]();
    int *coo_col = new int[nnz]();

    open_mtx_file(file, coo_col, coo_row, n, nnz);

    int *csr_row_ptr = new int[n + 1]();
    int *csr_col_ind = new int[nnz]();

    // std::cout << "\nCOO row:\t";
    // prt::arr(coo_row, nnz);
    // std::cout << "COO col:\t";
    // prt::arr(coo_col, nnz);

    coo2csr(csr_row_ptr, csr_col_ind, coo_row, coo_col, nnz, n, 0);

    delete[] coo_row;
    delete[] coo_col;

    std::cout << "\nCSR row_ptr:";
    prt::arr(csr_row_ptr, n + 1);
    std::cout << "CSR col_ind:";
    prt::arr(csr_col_ind, nnz);

    /* -------------------------------------------------------------------------- */


    /* ------------------------------- free memory ------------------------------ */

    delete[] csr_row_ptr;
    delete[] csr_col_ind;

    return 0;

}