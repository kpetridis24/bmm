/* -------------------------------------------------------------------------- */
/*                        parallel-masked-block-bmm.cpp                       */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <omp.h>

void parallelMaskedBlockBmm(int matIndF, int matIndA, int matIndB, int b)
{
    struct timeval timer;
    double t = -1;

    /* ------------------------------ read matrices ----------------------------- */

    int nF;
    int nA;
    int nB;
    int nnzF;
    int nnzA;
    int nnzB;
    csr F;
    csr A;
    csc B;

    timer = util::tic();

    read2csr(matIndF, nF, nnzF, F);
    read2csr(matIndA, nA, nnzA, A);
    read2csc(matIndB, nB, nnzB, B);

    t = util::toc(timer);
    std::cout << "\nReading of F, A and B completed\n" << "Reading time = " << t << " seconds" << std::endl;


    /* -------------------------------- blocking -------------------------------- */

    timer = util::tic();

    bcsr bcsrF;
    bcsr bcsrA; 
    bcsc bcscB;

    csr2bcsr(F, bcsrF, b);
    csr2bcsr(A, bcsrA, b);
    csc2bcsc(B, bcscB, b);

    t = util::toc(timer);
    std::cout << "\nBlocking of F, A and B completed\n" << "Blocking time = " << t << " seconds" << std::endl;

    util::delCsr(F);
    util::delCsr(A);
    util::delCsc(B);

    /* ----------------------------- block bmm test ----------------------------- */

    timer = util::tic();

    std::multimap<int, int> C;
    parallelMaskedBlockBmm(bcsrF, bcsrA, bcscB, C);

    t = util::toc(timer);
    std::cout << "\nParallel block-BMM completed\n" << "Parallel block-BMM time = " << t << " seconds" << std::endl;

    util::delBcsr(bcsrF);
    util::delBcsr(bcsrA);
    util::delBcsc(bcscB);

    std::vector<std::pair<int, int>> vecC;

    for (const auto& x : C) {
        vecC.push_back(std::pair<int, int> (x.first, x.second));
    }
    std::sort(vecC.begin(), vecC.end());

    /* ------------------------------ check result ------------------------------ */

    if (util::checkRes("C.mtx", vecC)) {
        std::cout << "\nTest passed\n";
    }
    else {
        std::cout << "\nTest failed\n";
    }
}

void parallelMaskedBlockBmm(bcsr &F, bcsr &A, bcsc &B, std::multimap <int, int> &C)
// masked boolean matrix multiplication F.*(A*B) using blocks
{
    if (A.n != B.m || A.m != F.m || B.n != F.n) {
        std::cout << "Dimensions error\n";
        exit(1);
    }

    if (A.b != B.b || A.b != F.b) {
        std::cout << "Block size error\n";
        exit(1);
    }

    int numBlockRowsF = F.m / F.b;

    for (int blockRowF = 0; blockRowF < numBlockRowsF; blockRowF++) {

        int _numOfNzBlocks = F.HL_bRowPtr[blockRowF + 1] - F.HL_bRowPtr[blockRowF];
        std::vector<std::multimap <int, int>> bRowC(_numOfNzBlocks);
        int startInd = F.HL_bRowPtr[blockRowF];

        #pragma omp parallel for schedule(dynamic)
        for (int indF = F.HL_bRowPtr[blockRowF]; indF < F.HL_bRowPtr[blockRowF + 1]; indF++) {

            int blockColF = F.HL_bColInd[indF];
            maskedBlockRowColMult(blockRowF, blockColF, F, A, B, bRowC[indF - startInd]); // add block to the block-row vector
        }

        for (int i = 0; i < _numOfNzBlocks; i++) {
            // the result block-rows of C have to be added in C sequentially to avoid data races
            int blockColF = F.HL_bColInd[F.HL_bRowPtr[blockRowF] + i];
            util::addCooBlockToMatrix(C, blockRowF, blockColF, A.b, bRowC[i]);
        }
    }
}