/* -------------------------------------------------------------------------- */
/*                                 headers.hpp                                */
/* -------------------------------------------------------------------------- */

#ifndef __HEADERS_HPP__
#define __HEADERS_HPP__

#include <iostream>
#include <bits/stdc++.h>
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
  int nnz;
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
    void cscMat(csc &M);
    void cooMat(coo &M);
    void map(std::multimap <int, int> m);
};

/* ---------------------------------- utils --------------------------------- */

namespace util
{
  struct timeval tic();
  static double toc(struct timeval begin);
  void blockOffsets(int blockInd, int *nzBlockIndex, int *blockNnzCounter, int b, int &LL_row_ptr_offset, int &LL_col_ind_offset);
  void addCooBlockToMatrix(int *M, int blockRow, int blockCol, int b, int &sizeM, std::multimap<int, int> &_mapC);
  bool checkRes(std::string checkGraph, coo &C);
  void initCsr(csr &M, int n, int nnz);
  void initCsc(csc &M, int n, int nnz);
  void initCoo(coo &M, int n, int nnz);
  void delCsr(csr &M);
  void delCsc(csc &M);
  void delCoo(coo &M);
  void delBcsr(bcsr &M);
  void delBcsc(bcsc &M);
};

/* --------------------------- blocking functions --------------------------- */

ret csr2bcsr(csr &M, bcsr &blM);
ret csr2bcsc(csr &M, bcsc &blM);

/* -------------------------------- block-bmm ------------------------------- */

bool rowColMult( int rowA, int colB,
                 bcsr &A, bcsc &B,
                 int LL_rowPtrOffsetA,
                 int LL_colIndOffsetA,
                 int LL_colPtrOffsetB,
                 int LL_rowIndOffsetB );
bool *blockRowColMult(int blockRowA, int blockColB, bcsr &A, bcsc &B);
void bbm( bcsr &A,
          bcsc &B,
          bool *_C,
          int LL_rowPtrOffsetA,
          int LL_colIndOffsetA,
          int LL_colPtrOffsetB,
          int LL_rowIndOffsetB );
void blockBmm(bcsr &A, bcsc &B);

ret2 *maskedBlockRowColMult( int blockRowA, int blockColB, 
                            bcsr &F, bcsr &A, bcsc &B, 
                            std::multimap<int, int> &_mapC );
void maskedBbm( bcsr &F,
                bcsr &A,
                bcsc &B,
                int LL_rowPtrOffsetF,
                int LL_colIndOffsetF,
                int LL_rowPtrOffsetA,
                int LL_colIndOffsetA,
                int LL_colPtrOffsetB,
                int LL_rowIndOffsetB,
                std::multimap <int, int> &_mapC );
ret2 maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B);
ret2 parallelMaskedBlockBmm(bcsr &F, bcsr &A, bcsc &B);

/* -------------------------------------------------------------------------- */

#endif