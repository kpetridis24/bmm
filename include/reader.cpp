/* -------------------------------------------------------------------------- */
/*                                 reader.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>

#include <headers.hpp>

void read2coo(int graphId, int &n, int &nnz, int &b, coo &A, coo &B)
{
    std::string graph;

    switch(graphId) {
        case 0:
            graph = "s6.mtx";
            // b = 2;
            b = 3;
            break;
        case 1:
            graph = "s12.mtx";
            b = 2;
            // b = 3;
            // b = 4;
            // b = 6;
            break;
        case 2:
            graph = "com-Youtube.mtx";
            b = 226978;
            // b = 113489;
            break;
        case 3:
            graph = "belgium_osm.mtx";
            b = 62665;
            break;
        case 4:
            graph = "dblp-2010.mtx";
            b = 23299;
            // b = 14182;
            break;
        case 5:
            graph = "as-Skitter.mtx";
            b = 48469;
            // b = 17857;
            break;
        default:
            exit(1);
    }

    std::string file = "graphs/" + graph;

    readMtxValues(file, n, nnz);

    util::initCoo(A, n, nnz);
    util::initCoo(B, n, nnz);

    openMtxFile(file, A.col, A.row, A.n, A.nnz);
    openMtxFile(file, B.col, B.row, B.n, B.nnz);
}

void read2csr(int graphId, int &n, int &nnz, int &b, csr &A, csc &B)
{
    std::string graph;

    switch(graphId) {
        case 0:
            graph = "s6.mtx";
            // b = 2;
            b = 3;
            break;
        case 1:
            graph = "s12.mtx";
            b = 2;
            // b = 3;
            // b = 4;
            // b = 6;
            break;
        case 2:
            graph = "com-Youtube.mtx";
            b = 226978;
            // b = 113489;
            break;
        case 3:
            graph = "belgium_osm.mtx";
            b = 62665;
            break;
        case 4:
            graph = "dblp-2010.mtx";
            b = 23299;
            // b = 14182;
            break;
        case 5:
            graph = "as-Skitter.mtx";
            b = 48469;
            // b = 17857;
            break;
        default:
            exit(1);
    }

    std::string file = "graphs/" + graph;

    readMtxValues(file, n, nnz);

    coo M;
    util::initCoo(M, n, nnz);

    openMtxFile(file, M.col, M.row, M.n, M.nnz);

    util::initCsr(A, n, nnz);
    util::initCsc(B, n, nnz);

    coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
    coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);

    util::delCoo(M);

    // prt::csrMat(A);
    // prt::cscMat(B);
}

void readMtxValues(std::string f, int &n, int &nnz)
{
    // Open the file:
    std::ifstream fin(f);
    // Declare variables:
    int M, N;
    // Ignore headers and comments:
    while (fin.peek() == '%') fin.ignore(2048, '\n');
    // Read defining parameters:
    fin >> M >> N >> nnz;
    if (M == N) {
        n = M;
    }
    else {
        std::cout << "not square" << std::endl;
    }
    fin.close();
}

// read coo in 0-based format
void openMtxFile(std::string f, int *row, int *col, int &n, int &nnz)
{
    // Open the file:
    std::ifstream fin(f);
    // Declare variables:
    int M, N;
    // Ignore headers and comments:
    while (fin.peek() == '%') fin.ignore(2048, '\n');
    // Read defining parameters:
    fin >> M >> N >> nnz;

    for (int i = 0; i < nnz; i++){
        fin >> row[i] >> col[i];
        row[i]--;
        col[i]--;
        // i++;
        // row[i]=col[i-1];
        // col[i]=row[i-1];
    }
    fin.close();
}

/* -------------------------------------------------------------------------- */
/*                                 COO to CSR                                 */
/* -------------------------------------------------------------------------- */

int coo2csr(
  int       * const row,       
  int       * const col,     
  int const * const cooRow,
  int const * const cooCol,  
  int const         nnz,      
  int const         n,         
  int const         isOneBased 
)
{
  // ----- cannot assume that input is already 0!
    for(int l = 0; l < n + 1; l++) row[l] = 0;

  // ----- find the correct column sizes
    for(int l = 0; l < nnz; l++)
        row[cooRow[l] - isOneBased]++;

  // ----- cumulative sum
    for(int i = 0, cumsum = 0; i < n; i++) {
        int temp = row[i];
        row[i] = cumsum;
        cumsum += temp;
    }
    row[n] = nnz;

  // ----- copy the row indices to the correct place
    for(int l = 0; l < nnz; l++) {
        int row_l = cooRow[l] - isOneBased;

        int dst = row[row_l];
        col[dst] = cooCol[l] - isOneBased;

        row[row_l]++;
    }

  // ----- revert the column pointers
    for(int i = 0, last = 0; i < n; i++) {
        int temp = row[i];
        row[i] = last;
        last = temp;
    }

    return n;
}

/* -------------------------------------------------------------------------- */