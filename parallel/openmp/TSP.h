#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <omp.h>

/*
 * Function that checks if a certain chromosome is
 * already in the array of chromosomes
 */
bool check_chromosome_openmp(int chromosome, int *chromosomes, int size) {
    if (chromosome == 0)
        return true;
    for (int i = 1; i < size; i++) {
        if (chromosome == chromosomes[i])
            return true;
    }
    return false;
}

/*
 * Function that generates a series of random chromosomes
 * to serve as first generation
 */
void generate_random_chromosomes_openmp(int *chromosomes, cities *c, int starting_point) {
    chromosomes[0] = starting_point;
    for (int i = 1; i < c->size; i++) {
        int city = rand() % c->size;

        while (check_chromosome_openmp(city, chromosomes, i))
            city = rand() % c->size;
        chromosomes[i] = city;
    }
    chromosomes[c->size] = starting_point;
}

/*
 * Function that constructs an individual by choosing the
 * minimum distance from each city when that particular
 * city is visited by the travelling salesman
 */
void minimum_chromosome_road_openmp(int *chromosomes, int start, cities *c) {
    int min, i, j, pos, city = start;

    chromosomes[0] = start;
    for (i = 0; i < c->size - 1; i++) {
        min = INT_MAX;
        for (j = 0; j < c->size; j++) {
            if (j != start && j != city &&
                c->roads[city][j] < min &&
                !check_chromosome_openmp(j, chromosomes, i)) {
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
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_openmp(individual **generation, cities *c,
    int start, int population_size, int start_index, int end_index) {
    int cost, i;
    int *chromosomes;
	
	for (i = start_index; i < end_index; i++) {
        cost = 0;
        chromosomes = generation[i]->chromosomes;
		for (int j = 0; j < c->size; j++) {
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
	for (i = 0; i < count_best; i++) {
		if (i < start_index)
			continue;
		if (i >= end_index)
			break;
        memcpy(next_generation[i]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
	for (i = current_index; i < current_index + count_best; i++) {
		if (i < start_index)
			continue;
		if (i >= end_index)
			break;
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
    for (i = current_index; i < current_index + count_best; i++) {
		if (i < start_index)
			continue;
		if (i >= end_index)
			break;
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

    for (i = current_index; i < population_size; i++) {
		if (i < start_index)
			continue;
		if (i >= end_index)
			break;
        generate_random_chromosomes_openmp(next_generation[i]->chromosomes, c, start);
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
 
    omp_set_num_threads(no_threads);

    t1 = omp_get_wtime();

    /* Step 1: Creating initial population */
    individual **current_generation = malloc(population_size * sizeof(individual*));
    individual **next_generation = malloc(population_size * sizeof(individual*));
    individual **auxiliary;

	#pragma omp parallel private(thread_id, start, end, i)
	{
		thread_id = omp_get_thread_num();
		start = thread_id * (population_size / no_threads);
		end = min((thread_id + 1) * (population_size / no_threads), population_size);
		
		for (i = start; i < end; i++) {
        	current_generation[i] = malloc(sizeof(individual));
    		current_generation[i]->fitness = 0;
    		current_generation[i]->position = i;
			current_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        	next_generation[i] = malloc(sizeof(individual));
    		next_generation[i]->fitness = 0;
    		next_generation[i]->position = i;
        	next_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        	for (j = 0; j < c->size + 1; j++) {
            	next_generation[i]->chromosomes[j] = 0;
    		}
    	}
	}


	minimum_chromosome_road_openmp(current_generation[0]->chromosomes, starting_point, c);
	for (int i = 1; i < population_size; i++) {
        generate_random_chromosomes_openmp(current_generation[i]->chromosomes, c, starting_point);
    }

	/* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
	#pragma omp parallel private(thread_id, start, end, i)
	{
    	for (i = 0; i < generations_no; i++) {
			thread_id = omp_get_thread_num();
			start = thread_id * (population_size / no_threads);
			end = min((thread_id + 1) * (population_size / no_threads), population_size);

        	/* Step 2: Calculating fitness */
			compute_generation_fitness_openmp(current_generation, c, starting_point,
				population_size, start, end);

			#pragma omp barrier
			#pragma omp single
			{
				/* Step 3: Sort in order of fitnesses */
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
			#pragma omp single
			{
        		/* Step 5: Switch to new generation */
        		auxiliary = current_generation;
        		current_generation = next_generation;
        		next_generation = auxiliary;
			}
		}

		compute_generation_fitness_openmp(current_generation, c, starting_point,
			population_size, start, end);
	}

	qsort(current_generation, population_size,
        sizeof(individual*), compare_individuals);

    t2 = omp_get_wtime();

    printf("Total execution time = %lf\n", t2 - t1);

    print_result_individual(current_generation, c);
}
