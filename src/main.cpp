/* -------------------------------------------------------------------------- */
/*                                  main.cpp                                  */
/* -------------------------------------------------------------------------- */

#include <headers.hpp>
#include <blocking.cpp>
#include <block-bmm.cpp>
#include <masked-block-bmm.cpp>
#include <parallel-masked-block-bmm.cpp>
#include <distributed-block-bmm.cpp>
#include <utils.cpp>
#include <reader.cpp>

int main(int argc, char **argv)
{
/* -------------------------------------------------------------------------- */
/*                                s12.mtx -> 0                                */
/*                                 F.mtx -> 1                                 */
/*                                 Α.mtx -> 2                                 */
/*                                 B.mtx -> 3                                 */
/* -------------------------------------------------------------------------- */

  int matIndF = 1;
  int matIndA = 2;
  int matIndB = 3;

/* ------------------------------- sequential ------------------------------- */

  // maskedBlockBmm(matIndF, matIndA, matIndB, argc, argv);

/* -------------------------------- parallel -------------------------------- */

  // parallelMaskedBlockBmm(matIndF, matIndA, matIndB, argc, argv);

/* ------------------------------- distributed ------------------------------ */

  bool isParallel = false;
  distributedBlockBmm(matIndF, matIndA, matIndB, isParallel, argc, argv);

/* -------------------------------------------------------------------------- */

  return 0;
}