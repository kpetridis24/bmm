/* -------------------------------------------------------------------------- */
/*                                blocking.cpp                                */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

ret csr2blocks( int *rowPtr, 
                int *colInd, 
                int N, 
                int nnz, 
                int b, 
                int *LL_bRowPtr, 
                int *LL_bColInd ) 
{
    int numBlocks = (N/b)*(N/b), emptyBlocks = 0;
    int *blockNnzCounter = new int[numBlocks + 1]();
    
    for(int i = 0; i < N; i++) 
        for(int j = rowPtr[i]; j < rowPtr[i + 1]; j++) 
            blockNnzCounter[(i / b) * (N / b) + (colInd[j] / b) + 1]++;
    
    for(int i = 0; i < numBlocks; i++) if (blockNnzCounter[i] == 0) emptyBlocks++;

    int t = 0, t2 = 0, blkPtrSize = numBlocks - emptyBlocks + 1;
    int *nzBlockIndex = new int[numBlocks]; 
    int *nzBlockIndex2 = new int[blkPtrSize];

    for(int i = 1; i < numBlocks + 1; i++) {

        blockNnzCounter[i] += blockNnzCounter[i - 1];
        if (blockNnzCounter[i] == blockNnzCounter[i - 1]) continue;
        nzBlockIndex[i - 1] = t++;
        nzBlockIndex2[t2++] = i - 1;
    }
    
    int cnt = 0;
    int *elementCounter = new int[numBlocks]();
    int blockIdx, relativeRowIdx, relativeColInd, colIndOffset;

    for (int i = 0; i < N; i++) {
        
        if (cnt == b) cnt = 0;

        for (int j = rowPtr[i]; j < rowPtr[i + 1]; j++) {

            blockIdx = (i / b) * (N / b) + (colInd[j] / b);
            colIndOffset = blockNnzCounter[blockIdx];
            LL_bColInd[colIndOffset + elementCounter[blockIdx]] = colInd[j] % b;
            elementCounter[blockIdx]++;
            LL_bRowPtr[nzBlockIndex[blockIdx] * (b + 1) + cnt + 1]++; 
        }
        cnt++;
    }

    delete[] elementCounter;

    int cumsum = 0;
    for (int l = 0; l < blkPtrSize; l++) {
        for (int v = l*(b+1); v < l*(b+1)+(b+1); v++) {
            cumsum += LL_bRowPtr[v];
            LL_bRowPtr[v] = cumsum;
        }
        cumsum = 0;
    }  

    prt::arr(blockNnzCounter, numBlocks+1);     //Non zeros of each block, thus externalBlockRowPtr
    prt::arr(nzBlockIndex, numBlocks);       //Non zero block indices, can be transformed to BCSR with the offsets
    
/* -------------------------------------------------------------------------- */
/*                                    TODO                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*                      pop nz blocks of blockNnzCounter                      */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*                                Low-Level CSR                               */
/* -------------------------------------------------------------------------- */
    
    std::cout << "\nLow-Level CSR\n";
    std::cout << "LL-b_rowPtr:\t";
    prt::arr(LL_bRowPtr, blkPtrSize * (b + 1));   //Inside blkRowPtr
    std::cout << "LL-bColInd:\t";
    prt::arr(LL_bColInd, nnz);

/* -------------------------------------------------------------------------- */
/*                                    B-COO                                   */
/* -------------------------------------------------------------------------- */

    int num_of_block_rows = N / b;
    int blocks_per_row = num_of_block_rows;
    int nnzb =  blkPtrSize;

    int *b_rows = new int[nnzb];
    int *b_cols = new int[nnzb];

    for (int i = 0; i < nnzb; i++) {
        b_rows[i] = nzBlockIndex2[i] / blocks_per_row;
        b_cols[i] = nzBlockIndex2[i] % blocks_per_row;
    }

    delete[] nzBlockIndex2;

    // prt::arr(b_rows, nnzb);
    // prt::arr(b_cols, nnzb);

/* -------------------------------------------------------------------------- */
/*                                    B-CSR                                   */
/* -------------------------------------------------------------------------- */
    
    int *HL_bRowPtr = new int[num_of_block_rows + 1];
    int *HL_bColInd = new int[nnzb];

    coo2csr(HL_bRowPtr, HL_bColInd, b_rows, b_cols, nnzb, num_of_block_rows, 0);

    std::cout << "\nHigh-Level B-CSR\n";
    std::cout << "HL-b_rowPtr:\t";
    prt::arr(HL_bRowPtr, num_of_block_rows + 1);
    std::cout << "HL-b_col_ind:\t";
    prt::arr(HL_bColInd, nnzb);

/* -------------------------------------------------------------------------- */

    delete[] b_rows;
    delete[] b_cols;

    ret _ret = {HL_bRowPtr, HL_bColInd, nzBlockIndex, blockNnzCounter};

    return _ret;
}

/* -------------------------------------------------------------------------- */
/*                                    main                                    */
/* -------------------------------------------------------------------------- */

/* ------------------------------ blocking test ----------------------------- */

// int b = 2;
// int numBlocks = (n / b) * (n / b);
// int LL_bRowPtrSize = numBlocks * (b + 1);
// int blocksPerRow = n / b;

// int *nzBlockIndex;
// int *blockNnzCounter;

// // Low-Level CSR
// int *LL_bRowPtr = new int[LL_bRowPtrSize]();
// int *LL_bColInd = new int[nnz]();

// // High-Level B-CSR
// int *HL_bRowPtr;
// int *HL_bColInd;

// // blocking
// ret _ret = csr2blocks(A.rowPtr, A.colInd, n, nnz, b, LL_bRowPtr, LL_bColInd);

// HL_bRowPtr = _ret.ret1;
// HL_bColInd = _ret.ret2;
// nzBlockIndex = _ret.ret3;
// blockNnzCounter = _ret.ret4;

/* ---------------------------- triCounting test ---------------------------- */

// int trNum = bCsrTriCount( LL_bRowPtr, 
//                         LL_bColInd, 
//                         HL_bRowPtr, 
//                         HL_bColInd,
//                         nzBlockIndex,
//                         blockNnzCounter, 
//                         n, 
//                         b );

// std::cout << "Num of triangles: " << trNum << std::endl;

/* ------------------------------- free memory ------------------------------ */

// delete[] LL_bRowPtr;
// delete[] LL_bColInd;
// delete[] HL_bRowPtr;
// delete[] HL_bColInd;
// delete[] nzBlockIndex;
// delete[] blockNnzCounter;