/* -------------------------------------------------------------------------- */
/*                                 reader.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>

void read2coo(int graphId, int &n, int &nnz, coo &M)
{
    std::string graph;

    switch(graphId) {
        case 0:
            graph = "s12.mtx";
            break;
        case 1:
            graph = "F.mtx";
            break;
        case 2:
            graph = "A.mtx";
            break;
        case 3:
            graph = "B.mtx";
            break;
        default:
            exit(1);
    }

    std::string file = "mtx/in/" + graph;

    readMtxValues(file, n, nnz);
    util::initCoo(M, n, n, nnz);
    openMtxFile(file, M.col, M.row, M.n, M.nnz);
    M.m = M.n;
}

std::string read2csr(int graphId, int &n, int &nnz, csr &A)
{
    std::string graph;

    switch(graphId) {
        case 0:
            graph = "s12.mtx";
            break;
        case 1:
            graph = "F.mtx";
            break;
        case 2:
            graph = "A.mtx";
            break;
        case 3:
            graph = "B.mtx";
            break;
        default:
            exit(1);
    }

    std::string file = "mtx/in/" + graph;

    readMtxValues(file, n, nnz);

    coo M;
    util::initCoo(M, n, n, nnz);

    openMtxFile(file, M.col, M.row, M.n, M.nnz);

    util::initCsr(A, n, n, nnz);

    coo2csr(A.rowPtr, A.colInd, M.row, M.col, A.nnz, A.n, 0);
    A.m = A.n;

    util::delCoo(M);

    return graph;
}

std::string read2csc(int graphId, int &n, int &nnz, csc &B)
{
    std::string graph;

    switch(graphId) {
        case 0:
            graph = "s12.mtx";
            break;
        case 1:
            graph = "F.mtx";
            break;
        case 2:
            graph = "A.mtx";
            break;
        case 3:
            graph = "B.mtx";
            break;
        default:
            exit(1);
    }

    std::string file = "mtx/in/" + graph;

    readMtxValues(file, n, nnz);

    coo M;
    util::initCoo(M, n, n, nnz);

    openMtxFile(file, M.col, M.row, M.n, M.nnz);

    util::initCsc(B, n, n, nnz);

    coo2csr(B.colPtr, B.rowInd, M.col, M.row, B.nnz, B.n, 0);
    B.m = B.n;

    util::delCoo(M);

    return graph;
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
    }
    fin.close();
}

int coo2csr(
  int       * const row,       
  int       * const col,     
  int const * const cooRow,
  int const * const cooCol,  
  int const         nnz,      
  int const         m,         
  int const         isOneBased 
)
{
  // ----- cannot assume that input is already 0!
    for(int l = 0; l < m + 1; l++) row[l] = 0;

  // ----- find the correct column sizes
    for(int l = 0; l < nnz; l++)
        row[cooRow[l] - isOneBased]++;

  // ----- cumulative sum
    for(int i = 0, cumsum = 0; i < m; i++) {
        int temp = row[i];
        row[i] = cumsum;
        cumsum += temp;
    }
    row[m] = nnz;

  // ----- copy the row indices to the correct place
    for(int l = 0; l < nnz; l++) {
        int row_l = cooRow[l] - isOneBased;

        int dst = row[row_l];
        col[dst] = cooCol[l] - isOneBased;

        row[row_l]++;
    }

  // ----- revert the column pointers
    for(int i = 0, last = 0; i < m; i++) {
        int temp = row[i];
        row[i] = last;
        last = temp;
    }

    return m;
}