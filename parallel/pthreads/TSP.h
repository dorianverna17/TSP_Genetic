#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <omp.h>
#include <pthread.h>
#include <math.h>

typedef struct info_t {
	int index;
	int population_size;
	int generations_no;
	int starting_point;
	int thread_id;
	int no_threads;
	int square_length;
	cities *c;
	individual **current_generation;
	individual **next_generation;
	individual **prev_generation;
	pthread_barrier_t *barrier;
} info;

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_pthreads(individual **generation, cities *c,
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
void mutate_generation_pthreads(individual **current_generation,
    individual **next_generation, cities *c, int start, int gen_no, int population_size,
	int start_index, int end_index) {
    int i, j, aux;
    int current_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;
	for (i = max(start_index, 0); i < min(count_best, end_index); i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
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
        generate_random_chromosomes(next_generation[i]->chromosomes, c, start);
    }
}

// function that merges the intervals from mergesort
// it is similar to the one implemented in the laboratory
void merge_intervals_pthreads(individual **source, int start, int mid, int end, individual **destination) {
	int iA = start;
	int iB = mid;
	int i;

	individual *aux;

	for (i = start; i < end; i++) {
		// here I implement the merging
		// the conditions are taken from the compare function
		// that was given as parameter to the qsort function
		// in the skel received from the APD team
		if (end == iB || (iA < mid && source[iB]->fitness > source[iA]->fitness)) {
			// destination[i] = source[iA];
			memcpy(destination[i], source[iA], sizeof(individual));
			iA++;
		} else if (end == iB || (iA < mid && source[iB]->fitness == source[iA]->fitness)) {
			// destination[i] = source[iA];
			memcpy(destination[i], source[iA], sizeof(individual));
			iA++;
		} else {
			// destination[i] = source[iB];
			memcpy(destination[i], source[iB], sizeof(individual));
			iB++;
		}
	}
}


void mergesort_parallel_pthreads(int thread_id, int actual_length, int square_length,
	individual **vNew, individual **v, int P, pthread_barrier_t *barrier) {
	int start_index, end_index;
	int start_local, end_local;
	int width;

	// declaring an auxiliary vector for the swap between
	// the vectors
	individual **aux;
	
	// we advance with the width of the vectors
	// that we are merging
	pthread_barrier_wait(barrier);
	// here are the steps performed by the mergesort
	// I gradually increase the width of the intervals which are merged
	// this part is similar to the one at the laboratory
	for (width = 1; width < actual_length; width = 2 * width) {
		// here I compute the starting and ending indexes 
		start_index = thread_id * (double) square_length / P;
		end_index = (thread_id + 1) * (double) square_length / P;
		// the the local index is greater than the squared length of
		// the generation to be sorted, hten the end index (local)
		// takes the value of this squared length
		start_local = (start_index / (2 * width)) * (2 * width);
		if (square_length > (end_index / (2 * width)) * (2 * width)) {
			end_local = (end_index / (2 * width)) * (2 * width);
		} else {
			end_local = square_length;
		}

		// here I iterate over the indexes which are specific to the thread
		// with the id thread_id
		for (int i = start_local; i < end_local; i = i + 2 * width) {
			// I have 3 cases in which I call the function that merges the intervals
			if (i + 2 * width > actual_length && i + width > actual_length) {
				// This is the case in which the size of the two intervals that
				// need to be merged ar actually overflowing, the size of the vector
				// overall is smaller. The position of the mid is also greater than
				// the actual size. This is why I set the mid and the end to be the
				// actual size of the vector 
				merge_intervals_pthreads(v, i, actual_length, actual_length, vNew);
			} else if (i + 2 * width > actual_length && i + width <= actual_length) {
				// The second case is just like the case above, but only the size of the
				// second interval is overflowing, the index of the mid is ok
				merge_intervals_pthreads(v, i, i + width, actual_length, vNew);
			} else {
				// here is the default case that has been treated in the laboratory as well
				merge_intervals_pthreads(v, i, i + width, i + 2 * width, vNew);
			}
		}

		// these barrier wait calls are for assuring that
		// all threads have finished execution of a part of code
		// needed by all of them later
		pthread_barrier_wait(barrier);
 
		// here I interchange the vectors so that I
		// get the result in v
		aux = v;
		v = vNew;
		vNew = aux;

		pthread_barrier_wait(barrier);
	}
}

void run_genetic_algorithm(info *information) {
	int index = information->index;
	int population_size = information->population_size;
	int generations_no = information->generations_no;
	int starting_point = information->starting_point;
	int thread_id = information->thread_id;
	int no_threads = information->no_threads;
	int square_length = information->square_length;
	
	int start = thread_id * (double) population_size / no_threads;
	int end = min((thread_id + 1) * (double) population_size / no_threads, population_size);
	int i, j;
	
	cities *c = information->c;

	individual **current_generation = information->current_generation;
	individual **next_generation = information->next_generation;
	individual **prev_generation = information->prev_generation;
	individual **auxiliary;

	pthread_barrier_t *barrier = information->barrier;
		
	for (i = start; i < end; i++) {
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

	pthread_barrier_wait(barrier);
	if (thread_id == 0) {	
		minimum_chromosome_road(current_generation[0]->chromosomes, starting_point, c);
		for (int i = 1; i < population_size; i++) {
       		generate_random_chromosomes(current_generation[i]->chromosomes, c, starting_point);
		}
	}

	/* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (i = 0; i < generations_no; i++) {
       	/* Step 2: Calculating fitness */
		pthread_barrier_wait(barrier);
		compute_generation_fitness_pthreads(current_generation, c, starting_point,
			population_size, start, end);
		pthread_barrier_wait(barrier);

		if (thread_id == 0) {
			qsort(current_generation, population_size,
        		sizeof(individual*), compare_individuals);
			/* transferred this line from the mutation function */
			generate_random_numbers(current_generation, c->size, population_size);
		}

		pthread_barrier_wait(barrier);
        /* Step 4: Selecting the best genes and mutating */
        mutate_generation_pthreads(current_generation, next_generation, c,
        	starting_point, i, population_size, start, end);

		pthread_barrier_wait(barrier);
        /* Step 5: Switch to new generation */
        auxiliary = current_generation;
        current_generation = next_generation;
        next_generation = auxiliary;
	}

	pthread_barrier_wait(barrier);
	compute_generation_fitness_pthreads(current_generation, c, starting_point,
		population_size, start, end);
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_parallel_pthreads(cities *c, int starting_point,
    int generations_no, int population_size, int no_threads) {
	pthread_t threads[no_threads];
	double t1, t2;

    individual **current_generation = malloc(population_size * sizeof(individual*));
    individual **next_generation = malloc(population_size * sizeof(individual*));
	individual **prev_generation = (individual **) malloc(population_size * sizeof(individual*));

	int res_pow = 1, square_length;
	int count_pow = 0;
	while (res_pow < population_size) {
		count_pow++;
		res_pow = pow(2, count_pow);
	}
	square_length = res_pow;

	//create the barrier
	pthread_barrier_t barrier;
	int err_barrier = pthread_barrier_init(&barrier, NULL, no_threads);
	if (err_barrier) {
		printf("Eroare la initializarea barierei\n");
		exit(-1);
	}

	info information[no_threads];

	for (int i = 0; i < no_threads; i++) {
        information[i].index = i;
		information[i].thread_id = i;
		information[i].no_threads = no_threads;
		information[i].population_size = population_size;
		information[i].generations_no = generations_no;
		information[i].current_generation = current_generation;
		information[i].next_generation = next_generation;
		information[i].prev_generation = prev_generation;
		information[i].square_length = square_length;
		information[i].starting_point = starting_point;
		information[i].c = c;
		information[i].barrier = &barrier;
    }

	t1 = omp_get_wtime();

	int err;
	for (int i = 0; i < no_threads; i++) {
		err = pthread_create(&threads[i], NULL, (void *)run_genetic_algorithm, &information[i]);
        if (err) {
            printf("Error when creating the thread\n");
            exit(-1);
        }
	}

	/* joining the threads when the algorithm is completed */
    for (int i = 0; i < no_threads; i++) {
		int err = pthread_join(threads[i], NULL);

		if (err) {
	  		printf("Eroare la asteptarea thread-ului %d\n", i);
	  		exit(-1);
		}
  	}

	pthread_barrier_destroy(&barrier);

    qsort(current_generation, population_size,
        sizeof(individual*), compare_individuals);

	t2 = omp_get_wtime();

	printf("Total execution time = %lf\n", t2 - t1);

    print_result_individual(current_generation, c);
}
