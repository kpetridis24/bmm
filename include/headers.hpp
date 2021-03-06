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
  int m;
  int n;
  int nnz;
} coo;

/* -------------------------------- CSR type -------------------------------- */

typedef struct
{
  int *rowPtr;
  int *colInd;
  int m;
  int n;
  int nnz;
} csr;

/* -------------------------------- CSC type -------------------------------- */

typedef struct
{
  int *colPtr;
  int *rowInd;
  int m;
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
  int m;
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
  int m;
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
  int size1;
  int size2;
  int size3;
  int size4;
} ret;

/* ------------------ masked block row-col mult return type ----------------- */

typedef struct
{
    int *M;
    int sizeM;
} ret2;

/* ----------------------------- read functions ----------------------------- */

void read2coo(int graphId, int &n, int &nnz, coo &M);
std::string read2csr(int graphId, int &n, int &nnz, csr &A);
std::string read2csc(int graphId, int &n, int &nnz, csc &B);
void readMtxValues(std::string f, int &n, int &nnz);
void openMtxFile(std::string f, int *row, int *col, int &n, int &nnz);
int coo2csr(
  int       * const row,       
  int       * const col,     
  int const * const cooRow,
  int const * const cooCol,  
  int const         nnz,      
  int const         m,         
  int const         isOneBased 
);

/* ----------------------------- print functions ---------------------------- */

namespace prt
{
    void arr(int *arr, int len);
    void mat(int **mat, int rows, int cols);
    void vec(std::vector<std::pair<int, int>> vec);
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
  void addCooBlockToMatrix(std::multimap<int, int> &mapM, int blockRow, int blockCol, int b, std::multimap<int, int> &_mapM);
  void removeCooRowOffsets(coo &M, int offset);
  void addCooRowOffsets(std::vector<std::pair<int, int>> &vecCooM, int *rowsM, int *colsM, int offset);
  void computeChunks(int* chunkSizes, int* chunkOffsets, int numProcesses, int numBlockRows);
  void initCsr(csr &M, int m, int n, int nnz);
  void initCsc(csc &M, int m, int n, int nnz);
  void initCoo(coo &M, int m, int n, int nnz);
  void delCsr(csr &M);
  void delCsc(csc &M);
  void delCoo(coo &M);
  void delBcsr(bcsr &M);
  void delBcsc(bcsc &M);
  bool checkRes(std::string checkGraph, std::vector <std::pair <int, int>> &vecC);
};

/* --------------------------- blocking functions --------------------------- */

ret csr2bcsr(csr &M, bcsr &blM);
ret csc2bcsc(csc &M, bcsc &blM);
void csr2bcsr(csr &M, bcsr &bcsrM, int b);
void csc2bcsc(csc &M, bcsc &bcscM, int b);

/* -------------------------------- block-bmm ------------------------------- */

bool rowColMult( int rowA, int colB,
                 bcsr &A, bcsc &B,
                 int LL_rowPtrOffsetA,
                 int LL_colIndOffsetA,
                 int LL_colPtrOffsetB,
                 int LL_rowIndOffsetB );

void blockBmm(bcsr &A, bcsc &B, std::multimap<int, int> &C);
void blockRowColMult( int blockRowA, int blockColB, 
                      bcsr &A, bcsc &B, 
                      std::multimap<int, int> &_mapC );
void bbm( bcsr &A,
          bcsc &B,
          int LL_rowPtrOffsetA,
          int LL_colIndOffsetA,
          int LL_colPtrOffsetB,
          int LL_rowIndOffsetB,
          std::multimap <int, int> &_C );

void maskedBlockBmm(bcsr &F, bcsr &A, bcsc &B, std::multimap<int, int> &mapC);
void maskedBlockRowColMult( int blockRowA, int blockColB, 
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

void parallelMaskedBlockBmm(bcsr &F, bcsr &A, bcsc &B, std::multimap <int, int> &C);

void maskedBlockBmm(int matIndF, int matIndA, int matIndB, int b);
void parallelMaskedBlockBmm(int matIndF, int matIndA, int matIndB, int b);

/* ------------------------------ mpi functions ----------------------------- */

void distributedBlockBmm(int matIndF, int matIndA, int matIndB, bool parallel, int b, int argc, char **argv);
double distributeCooMatrix(int numProcesses, int rank, coo &M, coo &_M, int matInd, int &b);
double broadcastCooMatrix(int numProcesses, int rank, coo &M, int matInd, int &b);
void bmmResultGather( int numProcesses,
                      int rank, 
                      int selfSize,
                      int &totalSize, 
                      int *rowsC, 
                      int *colsC, 
                      int *bmmResultRows, 
                      int *bmmResultCols,
                      std::vector <std::pair <int, int>> &bmmResultVec );

/* -------------------------------------------------------------------------- */

#endif