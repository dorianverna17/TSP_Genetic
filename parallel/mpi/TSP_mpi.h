#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"
#include "../../utils/mpi/mpi_utils.h"
#include "mpi.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>

#define ROOT 0

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_mpi(individual *generation, cities *c,
    int start, int start_pos, int end_pos, int rank) {
    int cost;
    int *chromosomes;
    for (int i = start_pos; i < end_pos; i++) {
        cost = 0;
        chromosomes = generation[i].chromosomes;

        for (int j = 0; j < c->size; j++) {
            cost += c->roads[chromosomes[j]][chromosomes[j + 1]];
        }
        generation[i].fitness = cost;
    }
}

/*
 * Function that generates two random numbers for
 * mutation swapping
 */
void generate_random_numbers_mpi(individual *gen, int size, int population_size, int rank) {

    // only one worker generates the numbers and broadcasts the values
    for (int i = 0; i < population_size; i++) {
        if (rank == ROOT) {
            gen[i].random_pos1 = 1 + rand() % (size - 1);
            gen[i].random_pos2 = 1 + rand() % (size - 1);
    
            while (gen[i].random_pos1 == gen[i].random_pos2) {
                gen[i].random_pos2 = 1 + rand() % (size - 1);
            }

            gen[i].random_pos3 = 1 + rand() % (size - 1);

            gen[i].random_pos4 = 1 + rand() % (size - 1);
            while (gen[i].random_pos4 == gen[i].random_pos3) {
                gen[i].random_pos4 = 1 + rand() % (size - 1);
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&(gen[i].random_pos1), 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&(gen[i].random_pos2), 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&(gen[i].random_pos3), 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&(gen[i].random_pos4), 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
    }
}


/*
 * Function that populates the next generation with the
 * first 40% individuals of the current generation, mutates
 * the first 30% and then crossover the first 30%
 */
void mutate_generation_mpi(individual *current_generation,
    individual *next_generation, cities *c, int start, int start_index, int end_index, int population_size, int rank) {
    int i, j, aux;
    int current_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;
    for (i = start_index; i < min(count_best, end_index); i++) {
        memcpy(next_generation[i].chromosomes,
            current_generation[i].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].fitness = 0;
    }

    current_index = count_best;
    
    // let's mutate the rest of them
    for (i = max(start_index, current_index); i < min(current_index + count_best, end_index); i++) {
        memcpy(next_generation[i].chromosomes,
            current_generation[i - current_index].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].random_pos1 =
            current_generation[i - current_index].random_pos1;
        next_generation[i].random_pos2 =
            current_generation[i - current_index].random_pos2;
        aux = next_generation[i].chromosomes[next_generation[i].random_pos1];
        next_generation[i].chromosomes[next_generation[i].random_pos1] =
            next_generation[i].chromosomes[next_generation[i].random_pos2];
        next_generation[i].chromosomes[next_generation[i].random_pos2] = aux;
    }

    current_index += count_best;

    // let's mutate the rest of them
    for (i = max(current_index, start_index); i < min(current_index + count_best, end_index); i++) {
        memcpy(next_generation[i].chromosomes,
            current_generation[i - current_index].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].random_pos1 =
            current_generation[i - current_index].random_pos1;
        next_generation[i].random_pos2 =
            current_generation[i - current_index].random_pos2;
        next_generation[i].random_pos3 =
            current_generation[i - current_index].random_pos3;
        next_generation[i].random_pos4 =
            current_generation[i - current_index].random_pos4;
        aux = next_generation[i].chromosomes[next_generation[i].random_pos1];
        next_generation[i].chromosomes[next_generation[i].random_pos1] =
            next_generation[i].chromosomes[next_generation[i].random_pos2];
        next_generation[i].chromosomes[next_generation[i].random_pos2] = aux;

        aux = next_generation[i].chromosomes[next_generation[i].random_pos3];
        next_generation[i].chromosomes[next_generation[i].random_pos3] =
            next_generation[i].chromosomes[next_generation[i].random_pos4];
        next_generation[i].chromosomes[next_generation[i].random_pos4] = aux;
    }

    current_index += count_best;

    for (i = max(start_index, current_index); i < min(population_size, end_index); i++) {
		memcpy(next_generation[i].chromosomes,
            current_generation[i - current_index].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].random_pos1 =
            current_generation[i - current_index].random_pos1;
        next_generation[i].random_pos2 =
            current_generation[i - current_index].random_pos2;
        next_generation[i].random_pos3 =
            current_generation[i - current_index].random_pos3;
        next_generation[i].random_pos4 =
            current_generation[i - current_index].random_pos4;
        aux = next_generation[i].chromosomes[next_generation[i].random_pos3];
        next_generation[i].chromosomes[next_generation[i].random_pos3] =
            next_generation[i].chromosomes[next_generation[i].random_pos2];
        next_generation[i].chromosomes[next_generation[i].random_pos2] = aux;

        aux = next_generation[i].chromosomes[next_generation[i].random_pos1];
        next_generation[i].chromosomes[next_generation[i].random_pos1] =
            next_generation[i].chromosomes[next_generation[i].random_pos4];
        next_generation[i].chromosomes[next_generation[i].random_pos4] = aux;
    }
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_parallel_mpi(cities *c, int starting_point, int generations_no, int population_size) {

	int rank, proc, global_pos, fitness;
    double t1, t2;
	int *auxiliary_chromosomes;
    
    individual *current_generation = NULL;
	individual *next_generation = NULL;
	individual *auxiliary = NULL;
    individual *received_current_generation = NULL;
    individual *received_next_generation = NULL;
    individual *aux_chromosomes;

	MPI_Init(NULL, NULL);
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    t1 = MPI_Wtime();

    /* Creating initial population */
    current_generation = malloc(population_size * sizeof(individual));
    next_generation = malloc(population_size * sizeof(individual));
    if (rank == ROOT)
        received_current_generation = malloc(population_size * sizeof(individual));
	/* Creating MPI Data_type for individuals*/
    MPI_Datatype mpi_individual_type;
    define_MPI_individual_type(&mpi_individual_type);

	
    for (int i = 0; i < population_size; i++) {
        /*
        * each worker has these 2 arrays and works only with a chunk of elements 
        * but needs access to the others elements for mutating
        */ 
        current_generation[i].fitness = 0;
        current_generation[i].position = i;
        current_generation[i].chromosomes = malloc((c->size + 1) * sizeof(int));
        next_generation[i].fitness = 0;
        next_generation[i].position = i;
        next_generation[i].chromosomes = calloc((c->size + 1), sizeof(int));

        /* array used when gathering the individuals and the chromosomes in ROOT */
        if (rank == 0)
            received_current_generation[i].chromosomes = malloc((c->size + 1) * sizeof(int));
    }

    /* auxiliary array used when receiving an array of chromosomes */
    aux_chromosomes = malloc((c->size + 1) * sizeof(int));

    /*
    * Used when gathering from workers - sendcounts[i] is the number of elements processed
    * by worker i
    */
    int local_ind_count = population_size / proc;
    int rem = population_size % proc;
    int *sendcounts = (int*)malloc(proc * sizeof(int));
    int *displs = (int*)malloc(proc * sizeof(int));
    int displ_aux = 0;
    for (int i = 0; i < proc; i++) {
        sendcounts[i] = local_ind_count;
        if(rem > 0) {
            sendcounts[i]++;
            rem--;
        }
        displs[i] = displ_aux;
        displ_aux += sendcounts[i];
    }
    

    /*
     * Allocate memory for the chromosomes of each generation, initialize them and send them to workers
    */
   	int rank_count = 0;
    int local_pos = 0;
   	if (rank == 0) {
		for (int i = 0; i < population_size; i++) {
			if (i == 0) {
				minimum_chromosome_road_mpi(current_generation[i].chromosomes, starting_point, c);
                local_pos++;
			} else {
				if (local_pos >= sendcounts[rank_count]) {
					local_pos = 0;
					rank_count++;
				}
				local_pos++;
				generate_random_chromosomes_mpi(current_generation[i].chromosomes, c, starting_point);
				if (rank_count != 0) {
					MPI_Isend(current_generation[i].chromosomes, (c->size + 1), MPI_INT, rank_count, i, MPI_COMM_WORLD, &request);
				}
			}
		}
    }
    
	if (rank != 0) {
		for (int j = 0; j < sendcounts[rank]; j++) {
			MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
        	MPI_Wait(&request, &status);
        	memcpy(current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
    	}
    }

    /* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (int i = 0; i < generations_no; i++) {
		/* Calculating fitness */
        global_pos = global_pos_calc(rank, 0, sendcounts, proc);
        compute_generation_fitness_mpi(current_generation, c, starting_point, global_pos, global_pos + sendcounts[rank], rank);
        MPI_Barrier(MPI_COMM_WORLD);

		/* send fitness from workers to ROOT */
        for (int j = 0; j < sendcounts[rank]; j++) {
            global_pos = global_pos_calc(rank, j, sendcounts, proc);
            MPI_Isend(&(current_generation[global_pos].fitness), 1, MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
        }

        if (rank == ROOT) {
            for (int j = 0; j < population_size; j++) {
                MPI_Irecv(&fitness, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                current_generation[status.MPI_TAG].fitness = fitness;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        
        /* send chromosomes from workers to ROOT */
        for (int j = 0; j < sendcounts[rank]; j++) {
            global_pos = global_pos_calc(rank, j, sendcounts, proc);
            MPI_Isend(current_generation[global_pos].chromosomes, (c->size + 1), MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
        }

        if (rank == ROOT) {
            for (int j = 0; j < population_size; j++) {
                MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, &status);
                memcpy(current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == ROOT) {
            /* Sort in order of fitnesses */
            qsort(current_generation, population_size,
                sizeof(individual), compare_individuals_mpi);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        /* Broadcasting the array of individuals and each array of chromosomes from each individual*/
        MPI_Bcast(current_generation, population_size, mpi_individual_type, ROOT, MPI_COMM_WORLD);
        for (int j = 0; j < population_size; j++) {
            MPI_Bcast(current_generation[j].chromosomes, c->size + 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        }


        generate_random_numbers_mpi(current_generation, c->size, population_size, rank);

        global_pos = global_pos_calc(rank, 0, sendcounts, proc);
		
		/* Selecting the best genes and mutating */
        mutate_generation_mpi(current_generation, next_generation, c,
            starting_point, global_pos, global_pos + sendcounts[rank], population_size, rank);
       
        /* Switch to new generation */
        auxiliary = current_generation;
        current_generation = next_generation;
        next_generation = auxiliary;

        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(&(current_generation[global_pos]), sendcounts[rank], mpi_individual_type, received_current_generation, sendcounts, displs, mpi_individual_type, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    /* send chromosomes from workers to ROOT */
    for (int j = 0; j < sendcounts[rank]; j++) {
        global_pos = global_pos_calc(rank, j, sendcounts, proc);
        MPI_Isend(current_generation[global_pos].chromosomes, (c->size + 1), MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
    }

    /* receiving the chromosomes */
    if (rank == ROOT) {
        for (int j = 0; j < population_size; j++) {
            MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            memcpy(received_current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
        }
    }

    if (rank == ROOT)
        memcpy(current_generation, received_current_generation, population_size * sizeof(individual));
    if (rank == 0) {

        compute_generation_fitness_mpi(current_generation, c, starting_point, 0, population_size, rank);

        qsort(current_generation, population_size,
            sizeof(individual), compare_individuals_mpi);
        
        t2 = MPI_Wtime();
        printf("Total execution time = %lf\n", t2 - t1);
        print_result_individual_mpi(current_generation, c);
    }
    

	MPI_Finalize();
}
