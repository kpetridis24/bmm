/* -------------------------------------------------------------------------- */
/*                                 reader.cpp                                 */
/* -------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>

void read_mtx_values(std::string f, int &n, int &nnz)
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
void open_mtx_file(std::string f, int *row, int *col, int &n, int &nnz)
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
  int const * const row_coo,
  int const * const col_coo,  
  int const         nnz,      
  int const         n,         
  int const         isOneBased 
)
{
  // ----- cannot assume that input is already 0!
    for(int l = 0; l < n+1; l++) row[l] = 0;


  // ----- find the correct column sizes
    for(int l = 0; l < nnz; l++)
        row[row_coo[l] - isOneBased]++;

  // ----- cumulative sum
    for(int i = 0, cumsum = 0; i < n; i++) {
        int temp = row[i];
        row[i] = cumsum;
        cumsum += temp;
    }
    row[n] = nnz;
  // ----- copy the row indices to the correct place
    for(int l = 0; l < nnz; l++) {
        int row_l = row_coo[l] - isOneBased;

        int dst = row[row_l];
        col[dst] = col_coo[l] - isOneBased;

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