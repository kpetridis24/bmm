/* -------------------------------------------------------------------------- */
/*                                 headers.hpp                                */
/* -------------------------------------------------------------------------- */

#ifndef __HEADERS_HPP__
#define __HEADERS_HPP__

#include <iostream>
#include <sys/time.h>

/* -------------------------------- COO type -------------------------------- */

typedef struct
{
  int *row;
  int *col;
  int n;
  int nnz;
} coo;

/* -------------------------------- CSR type -------------------------------- */

typedef struct
{
  int *rowPtr;
  int *colInd;
  int n;
  int nnz;
} csr;

/* -------------------------------- CSC type -------------------------------- */

typedef struct
{
  int *colPtr;
  int *rowInd;
  int n;
  int nnz;
} csc;

/* ------------------------------- B-CSR type ------------------------------- */

typedef struct
{
  int *LL_bRowPtr;
  int *LL_bColInd;
  int *HL_bRowPtr;
  int *HL_bColInd;
  int *nzBlockIndex;
  int *blockNnzCounter;
  int n;
  int b;
} bcsr;

/* ------------------------------- B-CSC type ------------------------------- */

typedef struct
{
  int *LL_bColPtr;
  int *LL_bRowInd;
  int *HL_bColPtr;
  int *HL_bRowInd;
  int *nzBlockIndex;
  int *blockNnzCounter;
  int n;
  int b;
} bcsc;

/* ---------------------- blocking function return type --------------------- */

typedef struct
{
  int *ret1;
  int *ret2;
  int *ret3; 
  int *ret4;
} ret;

/* ------------------ masked block row-col mult return type ----------------- */

typedef struct
{
    int *M;
    int sizeM;
} ret2;

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

namespace prt
{
    void arr(int *arr, int len);
    void mat(int **mat, int rows, int cols);
    void csrMat(csr &M);
    void cscMat(csr &M);
    void cooMat(csr &M);
};

/* ---------------------------------- utils --------------------------------- */

namespace util
{
  struct timeval tic();
  static double toc(struct timeval begin);
  void blockOffsets(int blockInd, int *nzBlockIndex, int *blockNnzCounter, int b, int &LL_row_ptr_offset, int &LL_col_ind_offset);
  void addCooElement(int row, int col, int *M, int &sizeM);
  bool searchCooElement(int row, int col, int *M, int &sizeM);
  void addCooBlockToMatrix(int *M, int *_M, int blockRow, int blockCol, int b, int &sizeM, int _sizeM);
  void initCsr(csr &M, int n, int nnz);
  void initCsc(csr &M, int n, int nnz);
  void initCoo(coo &M, int n, int nnz);
  void delCsr(csr &M);
  void delCsc(csr &M);
  void delCoo(csr &M);
  void delBcsr(coo &M);
  void delBcsc(coo &M);
};

/* --------------------------- blocking functions --------------------------- */

ret csr2bcsr(csr &M, bcsr &blM);
ret csr2bcsc(csr &M, bcsc &blM);

/* ----------------------------------- bmm ---------------------------------- */

bool rowColMult(int rowA, int colB, csr &A, csc &B);
void bmm(csr &A, csr &B, coo &C);
void maskedBmm(csr &F, csr &A, csc &B, coo &C);

/* -------------------------------- block-bmm ------------------------------- */

bool *blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B);
void bbm( bcsr &A,
          bcsc &B,
          bool *_C,
          int LL_rowPtrOffsetA,
          int LL_colIndOffsetA,
          int LL_colPtrOffsetB,
          int LL_rowIndOffsetB );
void blockBmm(bcsr &A, bcsc &B);
ret2 maskedBlockRowColMult(int blockRowA, int blockColB, bcsr &F, bcsr &A, bcsc &B);
void maskedBbm( bcsr &F,
                bcsr &A,
                bcsc &B,
                int *_C,
                int &_sizeC,
                int LL_rowPtrOffsetF,
                int LL_colIndOffsetF,
                int LL_rowPtrOffsetA,
                int LL_colIndOffsetA,
                int LL_colPtrOffsetB,
                int LL_rowIndOffsetB );
ret2 maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B);

/* -------------------------------------------------------------------------- */

#endif