#include "sequential/naive/TSP.h"
#include "sequential/genetic/TSP.h"
#include "parallel/openmp/TSP.h"
#include "parallel/pthreads/TSP.h"
#include "parallel/mpi/TSP_mpi.h"

#include <string.h>

#define NUM_THREADS 4
#define START_CITY 0

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: ./main algorithm_type input_file");
        printf("algorithm_type: sequential_naive, sequential_genetic, Pthreads, OpenMP, hibrid\n");
        return 1;
    }

    cities *c = read_file(argv[2]);

    if (strcmp(argv[1], "sequential_naive") == 0) {
        TSP_sequential_naive(c, 0);
    } else if (strcmp(argv[1], "sequential_genetic") == 0) {
        /* 1000 generations, 10000 individuals per generation */
        TSP_sequential_genetic(c, START_CITY, 1000, 10000);
    } else if (strcmp(argv[1], "parallel_openmp") == 0) {
        /* 1000 generations, 10000 individuals per generation */
        TSP_parallel_openmp(c, START_CITY, 1000, 10000, NUM_THREADS);
    }  else if (strcmp(argv[1], "parallel_pthreads") == 0) {
        /* 1000 generations, 10000 individuals per generation */
        TSP_parallel_pthreads(c, START_CITY, 1000, 10000, NUM_THREADS);
    } else if (strcmp(argv[1], "parallel_mpi") == 0) {
        /* 1000 generations, 1000 individuals per generation */
        TSP_parallel_mpi(c, START_CITY, 1000, 10000);
    }
}
