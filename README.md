# Parallel & Distributed Boolean Matrix Multiplication Using Blocks
This project exhibits the implementation of a fast, high-performance *BMM* algorithm for sparse matrices. Blocking is introduced as a practical way of accelerating *BMM* and facilitating its parallelism, while *MPI* and *OpenMP* are used to distribute and parallelize the computations respectively.

## Build instructions

Compile and run *sequential* version: 

```bash
make
```

Compile and run *parallel* version: 

```bash
make openmp
```

Compile and run *distributed* version: 

```bash
make mpi
```

Compile and run *hybrid* version: 

```bash
make hybrid
```

\* An *OpenMP* and *MPI* compatible compiler is required. 

\*\* Input datasets can be generated with `generator.m` and have to be placed in `mtx/in` folder.

\*\*\* In order for the tester to work, the result matrix `C.mtx` from `generator.m` have to be placed in  `mtx/out` folder.
