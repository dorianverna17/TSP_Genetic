#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_openmp(individual **generation, cities *c,
    int start, int population_size, int start_index, int end_index) {
    int cost, i, j;
    int *chromosomes;
	
	#pragma omp parallel for schedule(auto) private(i, j, cost, chromosomes)
	for (i = start_index; i < end_index; i++) {
        cost = 0;
        chromosomes = generation[i]->chromosomes;
		for (j = 0; j < c->size; j++) {
            cost += c->roads[chromosomes[j]][chromosomes[j + 1]];
        }
		generation[i]->fitness = cost;
    }
}

/*
 * Function that populates the next generation with the
 * first 40% individuals of the current generation, mutates
 * the first 30% and then crossover the first 30%
 */
void mutate_generation_openmp(individual **current_generation,
    individual **next_generation, cities *c, int start, int gen_no, int population_size,
	int start_index, int end_index) {
    int i, j, aux;
    int current_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;
	// #pragma omp parallel for schedule(auto) private(i)
	for (i = max(start_index, 0); i < min(count_best, end_index); i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
	// #pragma omp parallel for schedule(auto) private(i)
	for (i = max(start_index, current_index); i < min(current_index + count_best, end_index); i++) {
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

    current_index = current_index + count_best;

    // let's mutate the rest of them
	// #pragma omp parallel for schedule(auto) private(i)
    for (i = max(start_index, current_index); i < min(current_index + count_best, end_index); i++) {
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

    current_index = current_index + count_best;

    for (i = max(start_index, current_index); i < min(population_size, end_index); i++) {
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
        aux = next_generation[i]->chromosomes[next_generation[i]->random_pos3];
        next_generation[i]->chromosomes[next_generation[i]->random_pos3] =
            next_generation[i]->chromosomes[next_generation[i]->random_pos2];
        next_generation[i]->chromosomes[next_generation[i]->random_pos2] = aux;

        aux = next_generation[i]->chromosomes[next_generation[i]->random_pos1];
        next_generation[i]->chromosomes[next_generation[i]->random_pos1] =
            next_generation[i]->chromosomes[next_generation[i]->random_pos4];
        next_generation[i]->chromosomes[next_generation[i]->random_pos4] = aux;
    }
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_parallel_openmp(cities *c, int starting_point,
    int generations_no, int population_size, int no_threads) {

    int i, j, thread_id;
    double t1, t2;
	int start, end;

	t1 = omp_get_wtime();

    /* Step 1: Creating initial population */
    individual **current_generation = (individual **) malloc(population_size * sizeof(individual*));
    individual **next_generation = (individual **) malloc(population_size * sizeof(individual*));
	individual **prev_generation = (individual **) malloc(population_size * sizeof(individual*));
    individual **auxiliary;

	for (i = 0; i < population_size; i++) {
        current_generation[i] = malloc(sizeof(individual));
    	current_generation[i]->fitness = 0;
		current_generation[i]->position = i;
		current_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

		prev_generation[i] = malloc(sizeof(individual));
    	prev_generation[i]->fitness = 0;
		prev_generation[i]->position = i;
		prev_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        next_generation[i] = malloc(sizeof(individual));
    	next_generation[i]->fitness = 0;
    	next_generation[i]->position = i;
        next_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        for (j = 0; j < c->size + 1; j++) {
        	next_generation[i]->chromosomes[j] = 0;
    	}
	}

	minimum_chromosome_road(current_generation[0]->chromosomes, starting_point, c);
	for (int i = 1; i < population_size; i++) {
    	generate_random_chromosomes(current_generation[i]->chromosomes, c, starting_point);
	}

	omp_set_num_threads(no_threads);
	
	#pragma omp parallel private(thread_id, start, end, i, auxiliary)
	{
		thread_id = omp_get_thread_num();
		start = thread_id * (double) population_size / no_threads;
		end = min((thread_id + 1) * (double) population_size / no_threads, population_size);

		/* 
    	 * Complete the other steps until we reach the desired
    	 * number of generations
    	 */
    	for (i = 0; i < generations_no; i++) {
        	/* Step 2: Calculating fitness */
			#pragma omp barrier
			compute_generation_fitness_openmp(current_generation, c, starting_point,
				population_size, start, end);

			#pragma omp barrier
			// #pragma omp barrier
			if (thread_id == 0) {
				qsort(current_generation, population_size,
        			sizeof(individual*), compare_individuals);
				/* transferred this line from the mutation function */
				generate_random_numbers(current_generation, c->size, population_size);
			}

			#pragma omp barrier
        	/* Step 4: Selecting the best genes and mutating */
        	mutate_generation_openmp(current_generation, next_generation, c,
        		starting_point, i, population_size, start, end);

			#pragma omp barrier
			if (thread_id == 0) {
        		/* Step 5: Switch to new generation */
        		auxiliary = current_generation;
        		current_generation = next_generation;
        		next_generation = auxiliary;
			}
		}
	}

	t2 = omp_get_wtime();

	compute_generation_fitness_openmp(current_generation, c, starting_point,
		population_size, 0, population_size);

	qsort(current_generation, population_size,
        sizeof(individual*), compare_individuals);

    printf("Total execution time for OpenMP = %lf\n", t2 - t1);

    print_result_individual(current_generation, c);
}