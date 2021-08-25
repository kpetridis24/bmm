/* -------------------------------------------------------------------------- */
/*                                blocking.cpp                                */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

ret csr2bcsr(csr &M, bcsr &blM) 
{
    int numBlocks = (M.n / blM.b) * (M.n / blM.b), emptyBlocks = 0;
    int *blockNnzCounter = new int[numBlocks + 1]();
    
    for(int i = 0; i < M.n; i++) 
        for(int j = M.rowPtr[i]; j < M.rowPtr[i + 1]; j++) 
            blockNnzCounter[(i / blM.b) * (M.n / blM.b) + (M.colInd[j] / blM.b) + 1]++;
    
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
    int blockIdx, colIndOffset;

    for (int i = 0; i < M.n; i++) {
        
        if (cnt == blM.b) cnt = 0;

        for (int j = M.rowPtr[i]; j < M.rowPtr[i + 1]; j++) {

            blockIdx = (i / blM.b) * (M.n / blM.b) + (M.colInd[j] / blM.b);
            colIndOffset = blockNnzCounter[blockIdx];
            blM.LL_bColInd[colIndOffset + elementCounter[blockIdx]] = M.colInd[j] % blM.b;
            elementCounter[blockIdx]++;
            blM.LL_bRowPtr[nzBlockIndex[blockIdx] * (blM.b + 1) + cnt + 1]++; 
        }
        cnt++;
    }

    delete[] elementCounter;

    int cumsum = 0;
    for (int l = 0; l < blkPtrSize; l++) {
        for (int v = l*(blM.b + 1); v < l * ( blM.b + 1) + (blM.b + 1); v++) {
            cumsum += blM.LL_bRowPtr[v];
            blM.LL_bRowPtr[v] = cumsum;
        }
        cumsum = 0;
    }  

    // std::cout << "\nblockNnzCounter:";
    // prt::arr(blockNnzCounter, numBlocks+1);     //Non zeros of each block, thus externalBlockRowPtr
    // std::cout << "nzBlockIndex:";
    // prt::arr(nzBlockIndex, numBlocks);       //Non zero block indices, can be transformed to BCSR with the offsets
    
/* -------------------------------------------------------------------------- */
/*                                    TODO                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*                      pop nz blocks of blockNnzCounter                      */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*                                Low-Level CSR                               */
/* -------------------------------------------------------------------------- */
    
    // std::cout << "\nLow-Level CSR\n";
    // std::cout << "LL-bRowPtr:\t";
    // prt::arr(blM.LL_bRowPtr, blkPtrSize * (blM.b + 1));   //Inside blkRowPtr
    // std::cout << "LL-bColInd:\t";
    // prt::arr(blM.LL_bColInd, M.nnz);

/* -------------------------------------------------------------------------- */
/*                                    B-COO                                   */
/* -------------------------------------------------------------------------- */

    int num_of_block_rows = M.n / blM.b;
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

    // std::cout << "\nHigh-Level B-CSR\n";
    // std::cout << "HL-b_rowPtr:\t";
    // prt::arr(HL_bRowPtr, num_of_block_rows + 1);
    // std::cout << "HL-b_colInd:\t";
    // prt::arr(HL_bColInd, nnzb);

/* -------------------------------------------------------------------------- */

    delete[] b_rows;
    delete[] b_cols;

    ret _ret = {HL_bRowPtr, HL_bColInd, nzBlockIndex, blockNnzCounter};

    return _ret;
}

ret csr2bcsc(csc &M, bcsc &blM) 
{
    int numBlocks = (M.n / blM.b) * (M.n / blM.b), emptyBlocks = 0;
    int *blockNnzCounter = new int[numBlocks + 1]();
    
    for(int i = 0; i < M.n; i++) 
        for(int j = M.colPtr[i]; j < M.colPtr[i + 1]; j++) 
            blockNnzCounter[(i / blM.b) * (M.n / blM.b) + (M.rowInd[j] / blM.b) + 1]++;
    
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
    int blockIdx, rowIndOffset;

    for (int i = 0; i < M.n; i++) {
        
        if (cnt == blM.b) cnt = 0;

        for (int j = M.colPtr[i]; j < M.colPtr[i + 1]; j++) {

            blockIdx = (i / blM.b) * (M.n / blM.b) + (M.rowInd[j] / blM.b);
            rowIndOffset = blockNnzCounter[blockIdx];
            blM.LL_bRowInd[rowIndOffset + elementCounter[blockIdx]] = M.rowInd[j] % blM.b;
            elementCounter[blockIdx]++;
            blM.LL_bColPtr[nzBlockIndex[blockIdx] * (blM.b + 1) + cnt + 1]++; 
        }
        cnt++;
    }

    delete[] elementCounter;

    int cumsum = 0;
    for (int l = 0; l < blkPtrSize; l++) {
        for (int v = l*(blM.b + 1); v < l * ( blM.b + 1) + (blM.b + 1); v++) {
            cumsum += blM.LL_bColPtr[v];
            blM.LL_bColPtr[v] = cumsum;
        }
        cumsum = 0;
    }  

    // std::cout << "\nblockNnzCounter:";
    // prt::arr(blockNnzCounter, numBlocks+1);     //Non zeros of each block, thus externalBlockRowPtr
    // std::cout << "nzBlockIndex:";
    // prt::arr(nzBlockIndex, numBlocks);       //Non zero block indices, can be transformed to BCSR with the offsets
    
/* -------------------------------------------------------------------------- */
/*                                    TODO                                    */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/*                      pop nz blocks of blockNnzCounter                      */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/*                                Low-Level CSC                               */
/* -------------------------------------------------------------------------- */
    
    // std::cout << "\nLow-Level CSC\n";
    // std::cout << "LL-bColPtr:\t";
    // prt::arr(blM.LL_bColPtr, blkPtrSize * (blM.b + 1));   //Inside blkRowPtr
    // std::cout << "LL-bRowInd:\t";
    // prt::arr(blM.LL_bRowInd, M.nnz);

/* -------------------------------------------------------------------------- */
/*                                    B-COO                                   */
/* -------------------------------------------------------------------------- */

    int num_of_block_rows = M.n / blM.b;
    int blocks_per_row = num_of_block_rows;
    int nnzb =  blkPtrSize;

    int *b_rows = new int[nnzb];
    int *b_cols = new int[nnzb];

    for (int i = 0; i < nnzb; i++) {
        b_cols[i] = nzBlockIndex2[i] / blocks_per_row;  // TODO check
        b_rows[i] = nzBlockIndex2[i] % blocks_per_row;
    }

    delete[] nzBlockIndex2;

    // prt::arr(b_rows, nnzb);
    // prt::arr(b_cols, nnzb);

/* -------------------------------------------------------------------------- */
/*                                    B-CSC                                   */
/* -------------------------------------------------------------------------- */
    
    int *HL_bColPtr = new int[num_of_block_rows + 1];
    int *HL_bRowInd = new int[nnzb];

    coo2csr(HL_bColPtr, HL_bRowInd, b_cols, b_rows, nnzb, num_of_block_rows, 0);

    // std::cout << "\nHigh-Level B-CSC\n";
    // std::cout << "HL-bColPtr:\t";
    // prt::arr(HL_bColPtr, num_of_block_rows + 1);
    // std::cout << "HL-bRowInd:\t";
    // prt::arr(HL_bRowInd, nnzb);

/* -------------------------------------------------------------------------- */

    delete[] b_rows;
    delete[] b_cols;

    ret _ret = {HL_bColPtr, HL_bRowInd, nzBlockIndex, blockNnzCounter};

    return _ret;
}