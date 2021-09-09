/* -------------------------------------------------------------------------- */
/*                                blocking.cpp                                */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <headers.hpp>

ret csr2bcsr(csr &M, bcsr &blM) 
{
    blM.m = M.m;
    blM.n = M.n;

    int numBlocksM = (blM.m / blM.b) * (blM.n / blM.b);
    int LL_bRowPtrSizeM = numBlocksM * (blM.b + 1);

    // init Low-Level CSR
    blM.LL_bRowPtr = new int[LL_bRowPtrSizeM]();
    blM.LL_bColInd = new int[M.nnz]();


    int numBlockRows = M.m / blM.b;
    int blocksPerRow = M.n / blM.b;

    int numBlocks = (M.m / blM.b) * (M.n / blM.b);

    int emptyBlocks = 0;
    int nnzb = 0;
    int *blockNnzCounter = new int[numBlocks + 1]();
    bool *isNotEmpty = new bool[numBlocks]();
    
    for (int i = 0; i < M.m; i++) 
        for (int j = M.rowPtr[i]; j < M.rowPtr[i + 1]; j++) 
            blockNnzCounter[(i / blM.b) * blocksPerRow + (M.colInd[j] / blM.b) + 1]++;

    for (int i = 1; i < numBlocks + 1; i++) {
        if (blockNnzCounter[i] == 0) {
            emptyBlocks++;
        }
    }

    for (int i = 0; i < M.m; i++) 
        for (int j = M.rowPtr[i]; j < M.rowPtr[i + 1]; j++) 
            isNotEmpty[(i / blM.b) * blocksPerRow + (M.colInd[j] / blM.b)] = true;

    for (int i = 0; i < numBlocks; i++) {
        if (isNotEmpty[i])
            nnzb++;
    }

    delete[] isNotEmpty;

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

    for (int i = 0; i < M.m; i++) {
        
        if (cnt == blM.b) cnt = 0;

        for (int j = M.rowPtr[i]; j < M.rowPtr[i + 1]; j++) {

            blockIdx = (i / blM.b) * blocksPerRow + (M.colInd[j] / blM.b);
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

/* -------------------------------------------------------------------------- */
/*                                    B-COO                                   */
/* -------------------------------------------------------------------------- */

    int *b_rows = new int[nnzb];
    int *b_cols = new int[nnzb];

    for (int i = 0; i < nnzb; i++) {
        b_rows[i] = nzBlockIndex2[i] / blocksPerRow;
        b_cols[i] = nzBlockIndex2[i] % blocksPerRow;
    }

    delete[] nzBlockIndex2;

/* -------------------------------------------------------------------------- */
/*                                    B-CSR                                   */
/* -------------------------------------------------------------------------- */
    
    int *HL_bRowPtr = new int[numBlockRows + 1];
    int *HL_bColInd = new int[nnzb];

    coo2csr(HL_bRowPtr, HL_bColInd, b_rows, b_cols, nnzb, numBlockRows, 0);

    delete[] b_rows;
    delete[] b_cols;

    blM.HL_bRowPtr = HL_bRowPtr;
    blM.HL_bColInd = HL_bColInd;
    blM.nzBlockIndex = nzBlockIndex;
    blM.blockNnzCounter = blockNnzCounter;

    ret _ret = {HL_bRowPtr, HL_bColInd, nzBlockIndex, blockNnzCounter, numBlockRows+1,
                nnzb, numBlocks, numBlocks+1};
    return _ret;
}

ret csc2bcsc(csc &M, bcsc &blM) 
{

    blM.m = M.m;
    blM.n = M.n;

    int numBlocksM = (blM.m / blM.b) * (blM.n / blM.b);
    int LL_bColPtrSizeM = numBlocksM * (blM.b + 1);

    // init Low-Level CSC
    blM.LL_bColPtr = new int[LL_bColPtrSizeM]();
    blM.LL_bRowInd = new int[M.nnz]();

    int numBlockRows = M.m / blM.b;
    int blocksPerRow = M.n / blM.b;

    int numBlocks = (M.m / blM.b) * (M.n / blM.b), emptyBlocks = 0;
    int *blockNnzCounter = new int[numBlocks + 1]();
    int nnzb = 0;
    bool *isNotEmpty = new bool[numBlocks]();

    for(int i = 0; i < M.n; i++) 
        for(int j = M.colPtr[i]; j < M.colPtr[i + 1]; j++) {
            blockNnzCounter[(i / blM.b) * numBlockRows + (M.rowInd[j] / blM.b) + 1]++;
        }
    
    for(int i = 0; i < numBlocks; i++) {
        if (blockNnzCounter[i] == 0) {
            emptyBlocks++;
        }
    }

    for (int i = 0; i < M.n; i++) 
        for (int j = M.colPtr[i]; j < M.colPtr[i + 1]; j++) 
            isNotEmpty[(i / blM.b) * numBlockRows + (M.rowInd[j] / blM.b)] = true;

    for (int i = 0; i < numBlocks; i++) {
        if (isNotEmpty[i])
            nnzb++;
    }

    delete[] isNotEmpty;

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

/* -------------------------------------------------------------------------- */
/*                                    B-COO                                   */
/* -------------------------------------------------------------------------- */

    int *b_rows = new int[nnzb];
    int *b_cols = new int[nnzb];

    for (int i = 0; i < nnzb; i++) {
        b_cols[i] = nzBlockIndex2[i] / blocksPerRow;  // TODO check
        b_rows[i] = nzBlockIndex2[i] % blocksPerRow;
    }

    delete[] nzBlockIndex2;

/* -------------------------------------------------------------------------- */
/*                                    B-CSC                                   */
/* -------------------------------------------------------------------------- */
    
    int *HL_bColPtr = new int[numBlockRows + 1];
    int *HL_bRowInd = new int[nnzb];

    coo2csr(HL_bColPtr, HL_bRowInd, b_cols, b_rows, nnzb, numBlockRows, 0);

/* -------------------------------------------------------------------------- */

    delete[] b_rows;
    delete[] b_cols;

    blM.HL_bColPtr = HL_bColPtr;
    blM.HL_bRowInd = HL_bRowInd;
    blM.nzBlockIndex = nzBlockIndex;
    blM.blockNnzCounter = blockNnzCounter;

    ret _ret = {HL_bColPtr, HL_bRowInd, nzBlockIndex, blockNnzCounter, numBlockRows+1,
                nnzb, numBlocks, numBlocks + 1};

    return _ret;
}

void csr2bcsr(csr &M, bcsr &bcsrM, int b)
{
    bcsrM.b = b;

    // blocking
    csr2bcsr(M, bcsrM);
}

void csc2bcsc(csc &M, bcsc &bcscM, int b)
{
    bcscM.b = b;

    // blocking
    csc2bcsc(M, bcscM);
}