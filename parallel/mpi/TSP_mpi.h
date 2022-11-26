#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"
#include "mpi.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#define ROOT 0

/* 
 * Function that returns the index of an element from the global array
 * given the index from the local array in a process and the rank of the process 
 */
int global_index(int process_rank, int local_index, int global_size) {

}

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_mpi(individual **generation, cities *c,
    int start, int population_size) {
    int cost;
    int *chromosomes;

    for (int i = 0; i < population_size; i++) {
        cost = 0;
        chromosomes = generation[i]->chromosomes;
        for (int j = 0; j < c->size; j++) {
            cost += c->roads[chromosomes[j]][chromosomes[j + 1]];
        }
        generation[i]->fitness = cost;
    }
}

/*
 * Compare function used for sorting the individuals
 * from a certain generation
 */
int compare_individuals_mpi(const void *a, const void *b) {
    individual *individual_a = *(individual **) a;
    individual *individual_b = *(individual **) b;

    return individual_a->fitness - individual_b->fitness;
}

/*
 * Function that checks if a certain chromosome is
 * already in the array of chromosomes
 */
bool check_chromosome_mpi(int chromosome, int *chromosomes, int size) {
    if (chromosome == 0)
        return true;
    for (int i = 1; i < size; i++) {
        if (chromosome == chromosomes[i])
            return true;
    }
    return false;
}

/*
 * Function that prints the best individual in the resulted generation
 */
void print_result_individual_mpi(individual **gen, cities *c) {
    int *journey = gen[0]->chromosomes;

    printf("%d ", gen[0]->fitness);
    for (int i = 0; i < c->size; i++) {
        printf("%d -> ", journey[i] + 1);
    }
    printf("%d\n", journey[c->size] + 1);
}

void print_individual_mpi(int *chromosomes, int size, int fitness) {
    printf("%d ", fitness);

    for (int i = 0; i < size - 1; i++) {
        printf("%d -> ", chromosomes[i]);
    }
    printf("%d\n", chromosomes[size - 1]);
}

void print_generation_mpi(individual **gen, cities *c) {
    for (int i = 0; i < c->size; i++) {
        print_individual_mpi(gen[i]->chromosomes, c->size + 1, gen[i]->fitness);
    }
    printf("\n");
}

void generate_random_chromosomes_mpi(int *chromosomes, cities *c, int starting_point) {
    chromosomes[0] = starting_point;
    for (int i = 1; i < c->size; i++) {
        int city = rand() % c->size;

        while (check_chromosome_mpi(city, chromosomes, i))
            city = rand() % c->size;
        chromosomes[i] = city;
    }
    chromosomes[c->size] = starting_point;
}

/*
 * Function that generates two random numbers for
 * mutation swapping
 */
void generate_random_numbers_mpi(individual **gen, int size, int population_size) {
    for (int i = 0; i < population_size; i++) {
        gen[i]->random_pos1 = 1 + rand() % (size - 1);
        gen[i]->random_pos2 = 1 + rand() % (size - 1);
    
        while (gen[i]->random_pos1 == gen[i]->random_pos2) {
            gen[i]->random_pos2 = 1 + rand() % (size - 1);
        }

        gen[i]->random_pos3 = 1 + rand() % (size - 1);

        gen[i]->random_pos4 = 1 + rand() % (size - 1);
        while (gen[i]->random_pos4 == gen[i]->random_pos3) {
            gen[i]->random_pos4 = 1 + rand() % (size - 1);
        }
    }
}

/*
 * Function that constructs an individual by choosing the
 * minimum distance from each city when that particular
 * city is visited by the travelling salesman
 */
void minimum_chromosome_road_mpi(int *chromosomes, int start, cities *c) {
    int min, i, j, pos, city = start;

    chromosomes[0] = start;
    for (i = 0; i < c->size - 1; i++) {
        min = INT_MAX;
        for (j = 0; j < c->size; j++) {
            if (j != start && j != city &&
                c->roads[city][j] < min &&
                !check_chromosome_mpi(j, chromosomes, i)) {
                    min = c->roads[city][j];
                    pos = j;
            }
        }
        chromosomes[i + 1] = pos;
        city = pos;
    }
    chromosomes[i + 1] = start;
}

/*
 * Function that populates the next generation with the
 * first 40% individuals of the current generation, mutates
 * the first 30% and then crossover the first 30%
 */
void mutate_generation_mpi(individual **current_generation,
    individual **next_generation, cities *c, int start, int gen_no, int population_size) {
    int i, j, aux;
    int current_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;
    for (i = 0; i < count_best; i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
    generate_random_numbers_mpi(current_generation, c->size, population_size);
    for (i = current_index; i < current_index + count_best; i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i - current_index]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->random_pos1 =
            current_generation[i - current_index]->random_pos1;
        next_generation[i]->random_pos2 =
            current_generation[i - current_index]->random_pos2;
        aux = next_generation[i]->chromosomes[next_generation[i]->random_pos1];
        next_generation[i]->chromosomes[next_generation[i]->random_pos1] =
            next_generation[i]->chromosomes[next_generation[i]->random_pos2];
        next_generation[i]->chromosomes[next_generation[i]->random_pos2] = aux;
    }

    current_index = i;

    // let's mutate the rest of them
    for (i = current_index; i < current_index + count_best; i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i - current_index]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->random_pos1 =
            current_generation[i - current_index]->random_pos1;
        next_generation[i]->random_pos2 =
            current_generation[i - current_index]->random_pos2;
        next_generation[i]->random_pos3 =
            current_generation[i - current_index]->random_pos3;
        next_generation[i]->random_pos4 =
            current_generation[i - current_index]->random_pos4;
        aux = next_generation[i]->chromosomes[next_generation[i]->random_pos1];
        next_generation[i]->chromosomes[next_generation[i]->random_pos1] =
            next_generation[i]->chromosomes[next_generation[i]->random_pos2];
        next_generation[i]->chromosomes[next_generation[i]->random_pos2] = aux;

        aux = next_generation[i]->chromosomes[next_generation[i]->random_pos3];
        next_generation[i]->chromosomes[next_generation[i]->random_pos3] =
            next_generation[i]->chromosomes[next_generation[i]->random_pos4];
        next_generation[i]->chromosomes[next_generation[i]->random_pos4] = aux;
    }

    current_index = i;

    for (i = current_index; i < population_size; i++) {
        generate_random_chromosomes_mpi(next_generation[i]->chromosomes, c, start);
    }
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_parallel_mpi(cities *c, int starting_point, int generations_no, int population_size) {

	int rank, proc;
    double t1, t2; 

    individual aux_ind;
    MPI_Datatype individual_mpi, old_types[7];
    
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    t1 = MPI_Wtime();

    /* Step 1: Creating initial population */
	if (rank == ROOT) {
		individual **current_generation = malloc(population_size * sizeof(individual*));
		individual **next_generation = malloc(population_size * sizeof(individual*));
		individual **auxiliary;
	}
    
    individual **receive_generation;
	// needed to have contigous memory for the array

	individual *current_gen_aux = (individual *)malloc(population_size * population_size * sizeof(individual));
	if (rank == ROOT) {
		for (int i = 0; i < population_size; i++) {
				current_generation[i] = &(current_gen_aux[i * population_size]);
				current_generation[i]->fitness = 0;
				current_generation[i]->position = i;
				// current_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

				// if (i != 0) {
				// 	generate_random_chromosomes_mpi(current_generation[i]->chromosomes, c, starting_point);
				// } else {
				// 	minimum_chromosome_road_mpi(current_generation[i]->chromosomes, starting_point, c);
				// }

				next_generation[i] = malloc(sizeof(individual));
				next_generation[i]->fitness = 0;
				next_generation[i]->position = i;
				// next_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));
			
				// for (int j = 0; j < c->size + 1; j++) {
				// 	next_generation[i]->chromosomes[j] = 0;
				// }
			}
        
    }

    /* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (int i = 0; i < generations_no; i++) {
		MPI_Barrier(MPI_COMM_WORLD);

		//to do:
		// prepare scounts array for Scatterv + broadcast? chromosomes too
		MPI_Scatterv(current_generation, );
        /* Step 2: Calculating fitness */
        compute_generation_fitness_mpi(current_generation, c, starting_point, population_size);
        /* Step 3: Sort in order of fitnesses */
        qsort(current_generation, population_size,
            sizeof(individual*), compare_individuals_mpi);

        /* Step 4: Selecting the best genes and mutating */
        mutate_generation_mpi(current_generation, next_generation, c,
            starting_point, i, population_size);

        /* Step 5: Switch to new generation */
        auxiliary = current_generation;
        current_generation = next_generation;
        next_generation = auxiliary;
    }

    compute_generation_fitness_mpi(current_generation, c, starting_point, population_size);
    qsort(current_generation, population_size,
        sizeof(individual*), compare_individuals_mpi);
    
    t2 = MPI_Wtime();

    printf("Total execution time = %lf\n", t2 - t1);
    print_result_individual_mpi(current_generation, c);

	MPI_Finalize();
}
