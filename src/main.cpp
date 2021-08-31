/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <mpi.h>

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <block-bmm.cpp>
#include <masked-block-bmm.cpp>
// #include <parallel-masked-block-bmm.cpp>
#include <distributed-block-bmm.cpp>
#include <utils.cpp>
#include <reader.cpp>



int main(int argc, char **argv)
{
    /* -------------------------------------------------------------------------- */
    /*                        OpenMPI bmm distribution test                       */
    /* -------------------------------------------------------------------------- */

    /* --------------------------- Initialize OpenMPI --------------------------- */

    int numProcesses, pId;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &pId);
    // MPI_Datatype MPI_BCSC;
    MPI_Request req;
    MPI_Status stat;

    int arg1, arg2, arg3, arg4, arg5, arg6;
    int rxBuffer[6];
    bcsc blB;
    
    if(pId == 0) {
        
        struct timeval timer;
        double t = -1;

        /* ------------------------------- read matrix ------------------------------ */

        int n;
        int nnz;

        std::string graph = "s12.mtx";
        std::string file = "graphs/" + graph;

        readMtxValues(file, n, nnz);

        coo M;
        util::initCoo(M, n, nnz);

        openMtxFile(file, M.col, M.row, M.n, M.nnz);

        csr A;
        util::initCsr(A, n, nnz);
        csc B;
        util::initCsc(B, n, nnz);

        coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
        coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);

        util::delCoo(M);

        // prt::csrMat(A);
        // prt::cscMat(B);

        std::cout << "\nMatrix read successfully\nn = " << A.n << ", nnz = " << A.nnz << std::endl;

        /* ----------------------------------- s12 ---------------------------------- */

        // int b = 2;
        int b = 3;
        // int b = 4;
        // int b = 6;

        /* ------------------------------- com-Youtube ------------------------------ */

        // int b = 226978;
        // int b = 113489;
        
        /* -------------------------------- dblp-2010 ------------------------------- */

        // int b = 46598;
        // int b = 23299;
        // int b = 14182;
        // int b = 7091;
        // int b = 2026;
        // int b = 1013;

        /* ------------------------------- as-Skitter ------------------------------- */

        // int b = 242345;
        // int b = 89285;
        // int b = 48469;
        // int b = 17857;
        // int b = 12755;
        // int b = 2551;

        /* ------------------------------- belgium_osm ------------------------------ */

        // int b = 62665;

        /* --------------------------- bcsr blocking test --------------------------- */
        
        timer = util::tic();
        int numBlocks = (n / b) * (n / b);
        int LL_bRowPtrSize = numBlocks * (b + 1);

        bcsr blA;
        blA.n = A.n;
        blA.b = b;

        // init Low-Level CSR
        blA.LL_bRowPtr = new int[LL_bRowPtrSize]();
        blA.LL_bColInd = new int[nnz]();

        // blocking
        ret _ret = csr2bcsr(A, blA);

        blA.HL_bRowPtr = _ret.ret1;
        blA.HL_bColInd = _ret.ret2;
        blA.nzBlockIndex = _ret.ret3;
        blA.blockNnzCounter = _ret.ret4;

        t = util::toc(timer);
        std::cout << "\nBlocking A in B-CSR completed\n" << "Blocking time = " << t << " seconds" << std::endl;

        /* --------------------------- bcsc blocking test --------------------------- */

        timer = util::tic();

        int LL_bColPtrSize = numBlocks * (b + 1);

        // bcsc blB;
        blB.n = A.n;
        blB.b = b;

        // init Low-Level CSC
        blB.LL_bColPtr = new int[LL_bColPtrSize]();
        blB.LL_bRowInd = new int[nnz]();

        // blocking
        _ret = csr2bcsc(B, blB);

        blB.HL_bColPtr = _ret.ret1;
        blB.HL_bRowInd = _ret.ret2;
        blB.nzBlockIndex = _ret.ret3;
        blB.blockNnzCounter = _ret.ret4;

        /* ------------------- Send parameters to create bcsc type ------------------ */

        int sendBuffer[6] = {LL_bColPtrSize, nnz, _ret.size1, _ret.size2, _ret.size3, _ret.size4};

        arg1 = LL_bColPtrSize;
        arg2 = nnz;
        arg3 = _ret.size1;
        arg4 = _ret.size2;
        arg5 = _ret.size3;
        arg6 = _ret.size4;

        for(int p = 1; p < numProcesses; p++) {

            MPI_Isend(&sendBuffer, 6, MPI_INT, p, 0, MPI_COMM_WORLD, &req);
            MPI_Send(blB.LL_bColPtr, arg1, MPI_INT, p, 1, MPI_COMM_WORLD);
        }

        // std::cout<<blB.HL_bColPtr[1]<<std::endl;
        //for(int p = 1; p < numProcesses; p++) {
            
            // Send parameters to create bcsc mpi type
            // MPI_Isend(&sendBuffer, 6, MPI_INT, p, 0, MPI_COMM_WORLD, &req);
            // Send bcsc B array
            //MPI_Isend(&blB, 1, MPI_BCSC, p, 1, MPI_COMM_WORLD, &req);
        //}

        t = util::toc(timer);
        std::cout << "\nBlocking B in B-CSC completed\n" << "Blocking time = " << t << " seconds" << std::endl;
    }
    else {

        MPI_Recv(&rxBuffer, 6, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

        arg1 = rxBuffer[0];
        arg2 = rxBuffer[1];
        arg3 = rxBuffer[2];
        arg4 = rxBuffer[3];
        arg5 = rxBuffer[4];
        arg6 = rxBuffer[5];

        blB.LL_bColPtr = new int[arg1]();

        MPI_Recv(blB.LL_bColPtr, arg1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
        prt::arr(blB.LL_bColPtr, arg1);
    }

    








    // else {
    //     /* --------------------------- Receive parameters --------------------------- */
        
    //     bcsc blockB;
    //     MPI_Recv(&rxBuffer, 6, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

    //     arg1 = rxBuffer[0];
    //     arg2 = rxBuffer[1];
    //     arg3 = rxBuffer[2];
    //     arg4 = rxBuffer[3];
    //     arg5 = rxBuffer[4];
    //     arg6 = rxBuffer[5];
        
    //     // Create new MPI type
    //     //MPI_BCSC = createMpiTypeBcsc(arg1, arg2, arg3, arg4, arg5, arg6);
    //     const int nitems = 8;
    //     int blocklengths[8] = {arg1, arg2, arg3, arg4, arg5, arg6, 1, 1};
    //     MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    //     MPI_Datatype MPI_BCSC;
    //     MPI_Aint offsets[8];

    //     offsets[0] = offsetof(bcsc, LL_bColPtr);
    //     offsets[1] = offsetof(bcsc, LL_bRowInd);
    //     offsets[2] = offsetof(bcsc, HL_bColPtr);
    //     offsets[3] = offsetof(bcsc, HL_bRowInd);
    //     offsets[4] = offsetof(bcsc, nzBlockIndex);
    //     offsets[5] = offsetof(bcsc, blockNnzCounter);
    //     offsets[6] = offsetof(bcsc, n);
    //     offsets[7] = offsetof(bcsc, b);

    //     MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_BCSC);
    //     MPI_Type_commit(&MPI_BCSC);
    //     // Receive bcsc B array
    //     MPI_Recv(&blockB, 1, MPI_BCSC, 0, 1, MPI_COMM_WORLD, &stat);
    //     std::cout<<blockB.HL_bColPtr[1]<<std::endl;
    // }
    
    // MPI_BCSC = createMpiTypeBcsc(arg1, arg2, arg3, arg4, arg5, arg6);

    /* ------------------------ Creation of MPI_BCSC type ----------------------- */

    // const int nitems = 8;
    // int blocklengths[8] = {arg1, arg2, arg3, arg4, arg5, arg6, 1, 1};
    // MPI_Datatype types[8] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    // MPI_Datatype MPI_BCSC;
    // MPI_Aint offsets[8];

    // offsets[0] = offsetof(bcsc, LL_bColPtr);
    // offsets[1] = offsetof(bcsc, LL_bRowInd);
    // offsets[2] = offsetof(bcsc, HL_bColPtr);
    // offsets[3] = offsetof(bcsc, HL_bRowInd);
    // offsets[4] = offsetof(bcsc, nzBlockIndex);
    // offsets[5] = offsetof(bcsc, blockNnzCounter);
    // offsets[6] = offsetof(bcsc, n);
    // offsets[7] = offsetof(bcsc, b);

    // MPI_Type_create_struct(nitems, blocklengths, offsets, types, &MPI_BCSC);
    // MPI_Type_commit(&MPI_BCSC);

    // /* ---------------------------- Send BCSC array B --------------------------- */
    // //MPI_Bcast(&blB, 1, MPI_BCSC, 0, MPI_COMM_WORLD);
    
    // if(pId == 0) {

    //     for(int p = 0; p < numProcesses; p++)
    //         MPI_Send(&blB, 1, MPI_BCSC, p, 1, MPI_COMM_WORLD);
    // }
    // else {

    //     MPI_Recv(&blB, 1, MPI_BCSC, 0, 1, MPI_COMM_WORLD, &stat);
    //     //std::cout<<blB.HL_bColPtr[0]<<std::endl;
    // }
    //std::cout<<blB.HL_bColPtr[0]<<std::endl;
    

    // if(pId == 1) {

    //     bcsc blockB;
    //     MPI_Recv(&blockB, 1, MPI_BCSC, 0, 0, MPI_COMM_WORLD, &stat);
    // }
    
    // MPI_Bcast(&blB, 1, MPI_BCSC, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    /* ----------------------------- block bmm test ----------------------------- */

    // timer = util::tic();

    // std::multimap<int, int> C;

    // // blockBmm(blA, blB, C);
    // maskedBlockBmm(blA, blA, blB, C);
    // // ret2 ans = parallelMaskedBlockBmm(blA, blA, blB);

    // t = util::toc(timer);
    // std::cout << "\nBlock-BMM completed\n" << "Block-BMM time = " << t << " seconds" << std::endl;

    // std::vector<std::pair<int, int>> vecC;

    // for (const auto& x : C) {
    //   vecC.push_back(std::pair<int, int> (x.first, x.second));
    // }
    // std::sort(vecC.begin(), vecC.end());

    // // prt::vec(vecC);

    // /* ------------------------------ check result ------------------------------ */

    // if (util::checkRes(graph, vecC)) {
    //   std::cout << "\nTest passed\n";
    // }
    // else {
    //   std::cout << "\nTest failed\n";
    // }

    // /* ------------------------------- free memory ------------------------------ */

    // util::delCsr(A);
    // util::delCsc(B);
    // util::delBcsr(blA); 
    // util::delBcsc(blB);

  return 0;
}