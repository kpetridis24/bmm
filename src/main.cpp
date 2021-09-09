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

#include <headers.hpp>
#include <bmm.cpp>
#include <blocking.cpp>
#include <block-bmm.cpp>
#include <masked-block-bmm.cpp>
#include <parallel-masked-block-bmm.cpp>
#include <distributed-block-bmm.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main(int argc, char **argv)
{
  int matIndF = 6;
  int matIndA = 7;
  int matIndB = 8;

/* -------------------------------------------------------------------------- */
/*                                 sequential                                 */
/* -------------------------------------------------------------------------- */

  // maskedBlockBmm(matIndF, matIndA, matIndB, argc, argv);

/* -------------------------------------------------------------------------- */
/*                                  parallel                                  */
/* -------------------------------------------------------------------------- */

  // parallelMaskedBlockBmm(matIndF, matIndA, matIndB, argc, argv);

/* -------------------------------------------------------------------------- */
/*                                 distributed                                */
/* -------------------------------------------------------------------------- */

  bool isParallel = true;
  distributedBlockBmm(matIndF, matIndA, matIndB, isParallel, argc, argv);

  return 0;
}