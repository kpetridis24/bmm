/* -------------------------------------------------------------------------- */
/*                          distributed-block-bmm.cpp                         */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <mpi.h>

void distributedBlockBmm(int matIndF, int matIndA, int matIndB, bool isParallel, int b, int argc, char **argv)
{
    struct timeval timer;
    double tTotal = -1;
    double t1 = -1;
    double t2 = -1;
    double t3 = -1;
    double t4 = -1;
    double t5 = -1;
    double t6 = -1;

    int numProcesses, rank;

    /* ----------------- initialize MPI and declare COO matrices ---------------- */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    coo cooF;
    coo cooA;
    coo cooB;
    coo _cooF;
    coo _cooA;

    /* ---------------------- distribute A and broadcast B ---------------------- */

    t1 = distributeCooMatrix(numProcesses, rank, cooF, _cooF, matIndF, b);
    t2 = distributeCooMatrix(numProcesses, rank, cooA, _cooA, matIndA, b);
    t3 = broadcastCooMatrix(numProcesses, rank, cooB, matIndB, b);

    /* ------------------------ start timer in process 0 ------------------------ */

    if (rank == 0)
        timer = util::tic();

    /* ---------------- fix indices of matrix A - remove offsets ---------------- */

    int rowsPerChunk = _cooA.n / numProcesses;
    int chunkStartingRow = rowsPerChunk * rank;

    util::removeCooRowOffsets(_cooF, chunkStartingRow);
    util::removeCooRowOffsets(_cooA, chunkStartingRow);

    /* ------------------------ stop timer in process 0 ------------------------- */

    if (rank == 0)
        t4 = util::toc(timer);

    /* ------------------- convert F and A to CSR and B to CSC ------------------ */

    csr csrF;
    util::initCsr(csrF, _cooF.m, _cooF.n, _cooF.nnz);

    coo2csr(csrF.rowPtr, csrF.colInd, _cooF.row, _cooF.col, _cooF.nnz, _cooF.m, 0);

    csr csrA;
    util::initCsr(csrA, _cooA.m, _cooA.n, _cooA.nnz);

    coo2csr(csrA.rowPtr, csrA.colInd, _cooA.row, _cooA.col, _cooA.nnz, _cooA.m, 0);

    csc cscB;
    util::initCsc(cscB, cooB.m, cooB.n, cooB.nnz);

    coo2csr(cscB.colPtr, cscB.rowInd, cooB.col, cooB.row, cooB.nnz, cooB.n, 0);

    /* -------------------------- blocking of F Matrix -------------------------- */

    bcsr bcsrF;

    // blocking
    csr2bcsr(csrF, bcsrF, b);

    /* -------------------------- blocking of A Matrix -------------------------- */

    bcsr bcsrA;

    // blocking
    csr2bcsr(csrA, bcsrA, b);

    /* -------------------------- blocking of B Matrix -------------------------- */

    bcsc bcscB;

    // blocking
    csc2bcsc(cscB, bcscB, b);

    /* ------------------------ start timer in process 0 ------------------------ */

    if (rank == 0)
        timer = util::tic();

    /* -------------------------------- block BMM ------------------------------- */

    std::multimap <int, int> _cooC;

    if (!isParallel) {
        maskedBlockBmm(bcsrF, bcsrA, bcscB, _cooC);
    }
    else {
        parallelMaskedBlockBmm(bcsrF, bcsrA, bcscB, _cooC);
    }

    /* ------------------------ stop timer in process 0 ------------------------- */

    if (rank == 0)
        t5 = util::toc(timer);

    /* --------------------------- sort result matrix --------------------------- */

    std::vector <std::pair <int, int>> _vecCooC;

    for (auto &x : _cooC) {
      _vecCooC.push_back(std::pair <int, int> (x.first, x.second));
    }
    std::sort(_vecCooC.begin(), _vecCooC.end());

    /* ------------------------ start timer in process 0 ------------------------ */

    if (rank == 0)
        timer = util::tic();

    /* ------------------ fix indices of matrix C - add offsets ----------------- */

    int *_rowsC = new int[_vecCooC.size()];
    int *_colsC = new int[_vecCooC.size()];
    int *bmmResultRows;
    int *bmmResultCols;
    int selfSize = _vecCooC.size(), totalSize = 0;

    util::addCooRowOffsets(_vecCooC, _rowsC, _colsC, chunkStartingRow);
    std::vector <std::pair <int, int>> bmmResultVec;

    bmmResultGather(numProcesses, rank, selfSize, totalSize, _rowsC, _colsC, bmmResultRows,
                    bmmResultCols, bmmResultVec);

    /* ------------------------ stop timer in process 0 ------------------------- */

    if (rank == 0)
        t6 = util::toc(timer);

    // /* ------------------------------- free memory ------------------------------ */

    util::delCsr(csrF);
    util::delCsr(csrA);
    util::delCsc(cscB);
    util::delCoo(cooB);
    util::delCoo(_cooF);
    util::delCoo(_cooA);
    delete[] _rowsC;
    delete[] _colsC;
    delete[] bmmResultRows;
    delete[] bmmResultCols;

    /* ---------------------- compute total execution time ---------------------- */
    
    if (rank == 0) {
        tTotal = t1 + t2 + t3 + t4 + t5 + t6;
        std::cout << "\nDistributed block-BMM completed\n" << "Total time = " << tTotal << " seconds\n\n";
    }

    /* ------------------------------ MPI finalize ------------------------------ */

    MPI_Finalize();

    if (rank != 0)
        exit(0);

    /* ------------------------------ check result ------------------------------ */

    if (rank == 0) {
        if (util::checkRes("C.mtx", bmmResultVec)) {
            std::cout << "\nTest passed\n";
        }
        else {
            std::cout << "\nTest failed\n";
        }
    }
}

double distributeCooMatrix(int numProcesses, int rank, coo &M, coo &_M, int matInd, int &b)
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

    if (rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        read2coo(matInd, n, nnz, M);

        /* ------------------- compute num of block rows per chunk ------------------ */
        
        timer = util::tic();
        
        numBlockRows = n / b;

        if (numBlockRows % numProcesses) {
            std::cout << "numBlockRows % numProcesses != 0 \n";
            exit(1);
        }

        bRowsPerChunk = numBlockRows / numProcesses;
        rowsPerChunk = bRowsPerChunk * b;

        /* -------------------- compute num of elements per chunk ------------------- */

        for (int cnt = 0, elementCount = 0; cnt < M.nnz; cnt++, elementCount++) {

            currentElement = M.row[cnt];
            chunkInd = currentElement / rowsPerChunk;
            chunkSizes[chunkInd]++;
            chunkOffsets[chunkInd + 1] = chunkOffsets[chunkInd] + chunkSizes[chunkInd];
        }   

        selfChunkSize = chunkSizes[0];

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

    /* ---------------------------- distribute matrix --------------------------- */

    util::initCoo(_M, rowsPerChunk, n, selfChunkSize);

    MPI_Scatterv(M.row, chunkSizes, chunkOffsets, MPI_INT, _M.row,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);
    
    MPI_Scatterv(M.col, chunkSizes, chunkOffsets, MPI_INT, _M.col,
                 selfChunkSize, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        t = util::toc(timer);
        util::delCoo(M);
        return t;
    }

    delete[] chunkSizes;
    delete[] chunkOffsets;

    return 1;
}

double broadcastCooMatrix(int numProcesses, int rank, coo &M, int matInd, int &b)
{
    int n;
    int nnz;
    struct timeval timer;
    double t = -1;

    if(rank == 0) {

        /* ------------------------------- read matrix ------------------------------ */

        read2coo(matInd, n, nnz, M);

        /* ------------------------ start timer in process 0 ------------------------ */

        timer = util::tic();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank != 0) 
        util::initCoo(M, n, n, nnz);

    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.row, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(M.col, M.nnz, MPI_INT, 0, MPI_COMM_WORLD);

    t = util::toc(timer);

    if (rank == 0) {
        t = util::toc(timer);
        return t;
    }

    return 1;
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
    int *resultSizes = new int[numProcesses]();
    int *resultOffsets = new int[numProcesses]();
    int receivedSize, tag;

    if (rank != 0) {

        MPI_Send(&selfSize, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
    else {

        int cnt = 0;
        totalSize = selfSize;
        resultSizes[0] = selfSize;

        for (int p = 1; p < numProcesses; p++) {

            MPI_Recv(&receivedSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            tag = stat.MPI_TAG;
            resultSizes[tag] = receivedSize;
            totalSize += receivedSize;
        }

        for(int i = 1; i < numProcesses; i++) {
            resultOffsets[i] = resultOffsets[i - 1] + resultSizes[i - 1];
        }

        bmmResultRows = new int[totalSize];
        bmmResultCols = new int[totalSize];
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

    delete[] resultSizes;
    delete[] resultOffsets;
}