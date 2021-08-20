/* -------------------------------------------------------------------------- */
/*                                 headers.hpp                                */
/* -------------------------------------------------------------------------- */

#include <iostream>

/* ----------------------------- read functions ----------------------------- */

void read_mtx_values(std::string f, int &n, int &nnz);
void open_mtx_file(std::string f, int *row, int *col, int &n, int &nnz);
int coo2csr(
  int       * const row,       
  int       * const col,     
  int const * const row_coo,
  int const * const col_coo,  
  int const         nnz,      
  int const         n,         
  int const         isOneBased 
);

/* ----------------------------- print functions ---------------------------- */

namespace prt{
    void arr(int *arr, int len);
    void mat(int **mat, int rows, int cols);
};

/* ---------------------------------- utils --------------------------------- */

namespace util
{
    
};

/* -------------------------------------------------------------------------- */
