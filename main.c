#include "sequential/naive/TSP.h"
#include "sequential/genetic/TSP.h"

#include <string.h>

int main(int argc, char **argv) {
    if (argc < 3) {
        printf("Usage: ./main algorithm_type input_file");
        printf("algorithm_type: sequential_naive, sequential_genetic, Pthreads, OpenMP, hibrid");
        return 1;
    }

    cities *c = read_file(argv[2]);

    if (strcmp(argv[1], "sequential_naive") == 0) {
        TSP_sequential_naive(c, 0);
    } else if (strcmp(argv[1], "sequential_genetic") == 0) {
        TSP_sequential_genetic(c, 0, 1000000);
    }
}
