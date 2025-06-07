# Parallel Sequence Alignment

**Authors:** Ángela Garcinuño Feliciano and Lucas Massot

This project implements a parallel algorithm for sequence alignment. It includes:

- `simple.cpp`: A fully working parallel sequence alignment algorithm.
- `complex.cpp`: A partially working, more advanced implementation.

### Building the Project

To compile the code, run:

```bash
make
```

To execute the programs, run:

```bash
./<alg> <sequence_A_file> <sequence_B_file> <p>
```
With the inputs:

  - < alg > : Algorithm to use (simple / complex)
  
  - < sequence_A_file >: Name of the first sequence file (located in the sequences/ directory).

  - < sequence_B_file >: Name of the second sequence file (also in the sequences/ directory).

  - < p >: Number of processors to use.
