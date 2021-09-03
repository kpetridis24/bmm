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
    struct timeval timer;
    double t = -1;
    int matIndA = 0;
    int matIndB = 0;

  /* ----------------------- distributed block-BMM test ----------------------- */

    timer = util::tic();

    distributedBlockBmm(matIndA, matIndB, argc, argv);

    t = util::toc(timer);
    std::cout << "\nDistributed block-BMM completed\n" << "Total time = " << t << " seconds (pre-processing included)\n\n";

  return 0;
}