/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <headers.hpp>
#include <mpi.h>

void distributedBlockBmm(int matIndA, int matIndB, int argc, char **argv)
{
    struct timeval timer;
    double t = -1;

    int numProcesses, rank;
    int blockSizeA;
    int blockSizeB;

    /* ----------------- initialize MPI and declare COO matrices ---------------- */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    coo cooA;
    coo cooB;
    coo _cooA;

    /* ---------------------- distribute A and broadcast B ---------------------- */

    distributeCooMatrix(numProcesses, rank, cooA, _cooA, matIndA, blockSizeA);
    broadcastCooMatrix(numProcesses, rank, cooB, matIndB, blockSizeB);

    /* ---------------- fix indices of matrix A - remove offsets ---------------- */

    int rowsPerChunk = _cooA.n / numProcesses;
    int chunkStartingRow = rowsPerChunk * rank;

    util::removeCooRowOffsets(_cooA, chunkStartingRow);

    /* ---------------------- convert A to CSR and B to CSC --------------------- */

    csr csrA;
    util::initCsr(csrA, _cooA.m, _cooA.n, _cooA.nnz);

    coo2csr(csrA.rowPtr, csrA.colInd, _cooA.row, _cooA.col, _cooA.nnz, _cooA.m, 0);

    // prt::csrMat(csrA);

    csc cscB;
    util::initCsc(cscB, cooB.m, cooB.n, cooB.nnz);

    coo2csr(cscB.colPtr, cscB.rowInd, cooB.col, cooB.row, cooB.nnz, cooB.n, 0);

    // prt::cscMat(cscB);

    /* -------------------------- blocking of A Matrix -------------------------- */

    timer = util::tic();

    bcsr bcsrA;
    bcsrA.m = _cooA.m;
    bcsrA.n = _cooA.n;
    bcsrA.b = blockSizeA;

    int numBlocksA = (bcsrA.m / bcsrA.b) * (bcsrA.n / bcsrA.b);
    int LL_bRowPtrSize = numBlocksA * (bcsrA.b + 1);

    // init Low-Level CSR
    bcsrA.LL_bRowPtr = new int[LL_bRowPtrSize]();
    bcsrA.LL_bColInd = new int[_cooA.nnz]();

    // blocking
    ret _ret1 = csr2bcsr(csrA, bcsrA);

    bcsrA.HL_bRowPtr = _ret1.ret1;
    bcsrA.HL_bColInd = _ret1.ret2;
    bcsrA.nzBlockIndex = _ret1.ret3;
    bcsrA.blockNnzCounter = _ret1.ret4;
    int HL_bRowPtrSize = _ret1.size1;
    int HL_bColIndSize = _ret1.size2;
    int nzBlockIndexSizeA = _ret1.size3;
    int blockNnzCounterSizeA = _ret1.size4;


    t = util::toc(timer);
    // std::cout << "\nBlocking A in B-CSR completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    /* -------------------------- blocking of B Matrix -------------------------- */

    timer = util::tic();

    bcsc bcscB;
    bcscB.m = cooB.m;
    bcscB.n = cooB.n;
    bcscB.b = blockSizeB;

    int numBlocks = (bcscB.m / bcscB.b) * (bcscB.n / bcscB.b);
    int LL_bColPtrSize = numBlocks * (bcscB.b + 1);

    // init Low-Level CSC
    bcscB.LL_bColPtr = new int[LL_bColPtrSize]();
    bcscB.LL_bRowInd = new int[cooB.nnz]();

    // blocking
    ret _ret2 = csc2bcsc(cscB, bcscB);
    
    bcscB.HL_bColPtr = _ret2.ret1;
    bcscB.HL_bRowInd = _ret2.ret2;
    bcscB.nzBlockIndex = _ret2.ret3;
    bcscB.blockNnzCounter = _ret2.ret4;
    int HL_bColPtrSize = _ret2.size1;
    int HL_bRowIndSize = _ret2.size2;
    int nzBlockIndexSizeB = _ret2.size3;
    int blockNnzCounterSizeB = _ret2.size4;


    t = util::toc(timer);
    // std::cout << "\nBlocking B in B-CSC completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    /* -------------------------------- block BMM ------------------------------- */

    timer = util::tic();

    std::multimap <int, int> _cooC;
    // blockBmm(bcsrA, bcscB, _cooC);
    maskedBlockBmm(bcsrA, bcsrA, bcscB, _cooC);

    t = util::toc(timer);
    std::cout << "\nBlock-BMM completed\n" << "Block-BMM time = " << t << " seconds" << std::endl;

    std::vector <std::pair <int, int>> _vecCooC;

    timer = util::tic();

    for (auto &x : _cooC) {
      _vecCooC.push_back(std::pair <int, int> (x.first, x.second));
    }
    std::sort(_vecCooC.begin(), _vecCooC.end());

    t = util::toc(timer);
    // std::cout << "\nVector processing time = " << t << " seconds" << std::endl;

    /* ------------------ fix indices of matrix C - add offsets ----------------- */

    int *_rowsC = new int[_vecCooC.size()];
    int *_colsC = new int[_vecCooC.size()];
    int *bmmResultRows;
    int *bmmResultCols;
    int selfSize = _vecCooC.size(), totalSize = 0;

    util::addCooRowOffsets(_vecCooC, _rowsC, _colsC, chunkStartingRow);
    std::vector <std::pair <int, int>> bmmResultVec;
    
    // if (rank == 0) {
    //     std::cout << "Process " << rank << " result: \n";
    //     prt::arr(_rowsC, selfSize);
    //     prt::arr(_colsC, selfSize);
    // }

    bmmResultGather(numProcesses, rank, selfSize, totalSize, _rowsC, _colsC, bmmResultRows,
                    bmmResultCols, bmmResultVec);

    /* ------------------------------ MPI finalize ------------------------------ */

    MPI_Finalize();

    if (rank != 0)
        exit(0);

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);

    // /* ------------------------------ check result ------------------------------ */
    
    std::cout << "Distributed block-BMM result:\n";
    // prt::vec(bmmResultVec);

    // if (util::checkRes(matIndA, bmmResultVec)) {
    //   std::cout << "\nTest passed\n";
    // }
    // else {
    //   std::cout << "\nTest failed\n";
    // }
}

void distributeCooMatrix(int numProcesses, int rank, coo &M, coo &_M, int matInd, int &b)
{
    MPI_Request req;
    MPI_Status stat;

    int n;
    int nnz;
    int numBlockRows;
    int bRowsPerChunk;
    int rowsPerChunk; // = _m
    int currentElement;
    int selfChunkSize;
    int chunkInd;
    struct timeval timer;
    double t = -1;

    int *chunkSizes = new int[numProcesses]();
    int *chunkOffsets = new int[numProcesses]();

    /* -------------------------------- process 0 ------------------------------- */

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */
        
        timer = util::tic();

        read2coo(matInd, n, nnz, b, M);
        
        // std::cout << "\nMatrix read successfully\nn = " << M.n << ", nnz = " << M.nnz << std::endl;
        t = util::toc(timer);
        // std::cout << "\nReading time = " << t << std::endl;
        // prt::cooMat(M);

        /* ------------------- compute num of block rows per chunk ------------------ */
        
        timer = util::tic();
        
        numBlockRows = n / b;

        if (numBlockRows % numProcesses) {
            std::cout << "numBlockRows % numProcesses != 0 \n";
            exit(1);
        }

        bRowsPerChunk = numBlockRows / numProcesses;
        rowsPerChunk = bRowsPerChunk * b;

        // std::cout << "b = " << b << std::endl;
        // std::cout << "numProcesses = " << numProcesses << std::endl;
        // std::cout << "numBlockRows = " << numBlockRows << std::endl;
        // std::cout << "bRowsPerChunk = " << bRowsPerChunk << std::endl;

        /* -------------------- compute num of elements per chunk ------------------- */

        for (int cnt = 0, elementCount = 0; cnt < M.nnz; cnt++, elementCount++) {

            currentElement = M.row[cnt];
            chunkInd = currentElement / rowsPerChunk;
            // std::cout << "chunkId=" << chunkInd << std::endl;

            chunkSizes[chunkInd]++;
            chunkOffsets[chunkInd + 1] = chunkOffsets[chunkInd] + chunkSizes[chunkInd];
        }   

        selfChunkSize = chunkSizes[0];

        // prt::arr(chunkSizes, numProcesses);
        // prt::arr(chunkDisplacements, numProcesses);

        /* ---------- send array sizes to allocate memory on all processes ---------- */
        
        for (int p = 1; p < numProcesses; p++) 
            MPI_Send(&chunkSizes[p], 1, MPI_INT, p, 0, MPI_COMM_WORLD);
        
    }
    else {
        /* --------------------------- receive array sizes -------------------------- */

        MPI_Recv(&selfChunkSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

    }

    /* --------------------------- broadcast variables -------------------------- */

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowsPerChunk, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // std::cout << "chunk " << rank << "\t_m = " << rowsPerChunk << "\tn = " << n << "\tb = " << b << std::endl;

    /* ---------------------------- distribute matrix --------------------------- */

    util::initCoo(_M, rowsPerChunk, n, selfChunkSize);

    MPI_Scatterv(M.row, chunkSizes, chunkOffsets, MPI_INT, _M.row,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Scatterv(M.col, chunkSizes, chunkOffsets, MPI_INT, _M.col,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        t = util::toc(timer);
        // std::cout << "Distribution of matrix time = " << t << std::endl;
        util::delCoo(M);
    }

    // prt::cooMat(_M);
}

void broadcastCooMatrix(int numProcesses, int rank, coo &M, int matInd, int &b)
{
    int n;
    int nnz;
    struct timeval timer;
    double t = -1;

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        timer = util::tic();

        read2coo(matInd, n, nnz, b, M);
        
        // std::cout << "\nMatrix read successfully\nn = " << M.n << ", nnz = " << M.nnz << std::endl;
        t = util::toc(timer);
        // std::cout << "\nReading time = " << t << std::endl;
        // prt::cooMat(M);
    }

    timer = util::tic();

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank != 0) 
        util::initCoo(M, n, n, nnz);

    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.row, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.col, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);

    t = util::toc(timer);
    // std::cout << "Broadcast of matrix time = " << t << std::endl;

    // prt::cooMat(M);
}

void bmmResultGather( int numProcesses,
                      int rank, 
                      int selfSize,
                      int &totalSize, 
                      int *rowsC, 
                      int *colsC, 
                      int *bmmResultRows, 
                      int *bmmResultCols,
                      std::vector <std::pair <int, int>> &bmmResultVec )
{
    MPI_Status stat;
    int *resultSizes = new int[numProcesses];
    int *resultOffsets = new int[numProcesses];
    int receivedSize, tag;

    if (rank != 0) {

        MPI_Send(&selfSize, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
    else {

        int cnt = 0;
        totalSize = selfSize;
        resultSizes[0] = selfSize;
        resultOffsets[0] = 0;
        
        for (int p = 1; p < numProcesses; p++) {

            MPI_Recv(&receivedSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            tag = stat.MPI_TAG;
            resultSizes[tag] = receivedSize;
            resultOffsets[tag + 1] = resultOffsets[tag] + resultSizes[tag];
            totalSize += receivedSize;
        }

        bmmResultRows = new int[totalSize];
        bmmResultCols = new int[totalSize];
        // prt::arr(resultSizes, numProcesses);
        // prt::arr(resultOffsets, numProcesses);
    }

    MPI_Gatherv(rowsC, selfSize, MPI_INT, bmmResultRows, resultSizes, resultOffsets,
                MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Gatherv(colsC, selfSize, MPI_INT, bmmResultCols, resultSizes, resultOffsets,
                MPI_INT, 0, MPI_COMM_WORLD);

    // Transfer to vector in order to use tester
    if (rank == 0) {

        for (int i = 0; i < totalSize; i++) 
            bmmResultVec.push_back(std::pair <int, int> (bmmResultRows[i], bmmResultCols[i]));
        std::sort(bmmResultVec.begin(), bmmResultVec.end());
    } 
}