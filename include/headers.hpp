/* -------------------------------------------------------------------------- */
/*                                 headers.hpp                                */
/* -------------------------------------------------------------------------- */

#ifndef __HEADERS_HPP__
#define __HEADERS_HPP__

#include <iostream>

/* ---------------------- blocking function return type --------------------- */

typedef struct
{
    int *ret1;
    int *ret2;
    int *ret3;
    int *ret4;
} ret;

/* ----------------------------- read functions ----------------------------- */

void readMtxValues(std::string f, int &n, int &nnz);
void openMtxFile(std::string f, int *row, int *col, int &n, int &nnz);
int coo2csr(
  int       * const row,       
  int       * const col,     
  int const * const cooRow,
  int const * const cooCol,  
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

#endif