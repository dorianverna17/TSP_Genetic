#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"
#include "../../utils/mpi/mpi_utils.h"
#include "mpi.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

#define ROOT 0

typedef struct info_thr {
	int index;
	int population_size;
	int generations_no;
	int starting_point;
	int thread_id;
	int no_threads;
	cities *c;
	individual *current_generation;
	individual *next_generation;
	pthread_barrier_t *barrier;
} info_thr;

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_mpi_pthreads(individual *generation, cities *c,
    int start, int population_size, int no_threads) {
    int cost, start_index, end_index, i, thread_id, j;
    int *chromosomes;

    start_index = thread_id * (double) population_size / no_threads;
    end_index = min((thread_id + 1) * (double) population_size / no_threads, population_size);

    for (i = start_index; i < end_index; i++) {

        cost = 0;
        chromosomes = generation[i].chromosomes;

        for (j = 0; j < c->size; j++) {
            cost += c->roads[chromosomes[j]][chromosomes[j + 1]];
        }
        generation[i].fitness = cost;
    }
}


/*
 * Function that generates two random numbers for
 * mutation swapping
 */
void generate_random_numbers_mpi_pthreads(individual *gen, int size, int population_size) {

    // only one worker generates the numbers and broadcasts the values
    for (int i = 0; i < population_size; i++) {
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
}

void pthreads_mutate(info_thr *information) {
    int index = information->index;
	int population_size = information->population_size;
	int thread_id = information->thread_id;
	int no_threads = information->no_threads;
	
	int start = thread_id * (double) population_size / no_threads;
	int end = min((thread_id + 1) * (double) population_size / no_threads, population_size);
	
	cities *c = information->c;
    individual *current_generation = information->current_generation;
	individual *next_generation = information->next_generation;
    
    int current_index, aux, i, j;

    int count_best = (population_size * 3) / 10;

    for (i = start; i < min(count_best, end); i++) {
        memcpy(next_generation[i].chromosomes,
            current_generation[i].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
    for (i = max(start, current_index); i < min(current_index + count_best, end); i++) {
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

    current_index = i;

    // let's mutate the rest of them
    for (i = max(start, current_index); i < min(current_index + count_best, end); i++) {
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

    current_index = i;

    for (i = max(start, current_index); i < min(population_size, end); i++) {
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
void TSP_parallel_mpi_pthreads(cities *c, int starting_point, int generations_no, int population_size, int no_threads) {

	int rank, proc, global_pos;
    double t1, t2; 
    
    individual *current_generation = NULL;
	individual *next_generation = NULL;
	individual *auxiliary = NULL;
    individual *received_current_generation = NULL;
    int *aux_chromosomes;
    int provided;

	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    /* Create MPI Data_type for individuals*/
    MPI_Datatype mpi_individual_type;
    define_MPI_individual_type(&mpi_individual_type);

    /* pthreads initialization */
    pthread_t threads[no_threads];
	info_thr information[no_threads];
    for (int i = 0; i < no_threads; i++) {
        information[i].index = i;
		information[i].thread_id = i;
		information[i].no_threads = no_threads;
		information[i].population_size = population_size;
		information[i].generations_no = generations_no;
		information[i].starting_point = starting_point;
		information[i].c = c;
    }

    t1 = omp_get_wtime();

    /* Step 1: Creating initial population in the ROOT worker*/
    if (rank == ROOT) {
        current_generation = malloc(population_size * sizeof(individual));
        next_generation = malloc(population_size * sizeof(individual));
        for (int i = 0; i < population_size; i++) {
            current_generation[i].fitness = 0;
            current_generation[i].position = i;
            current_generation[i].chromosomes = malloc((c->size + 1) * sizeof(int));

            if (i == 0) {
                minimum_chromosome_road_mpi(current_generation[i].chromosomes, starting_point, c);
            } else {
                generate_random_chromosomes_mpi(current_generation[i].chromosomes, c, starting_point);
            }

            next_generation[i].fitness = 0;
            next_generation[i].position = i;
            next_generation[i].chromosomes = calloc((c->size + 1), sizeof(int));
        }
       
    }
    

    /*
    * Used when Scattering and Gathering from workers - sendcounts[i] is the number of elements processed
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
    
    /* the arrays all workers will use - each has only a chunk of elements from the original */
    received_current_generation = (individual*)malloc(sendcounts[rank] * sizeof(individual));

    /*
     * Allocate memory for the chromosomes in the workers
    */
    for (int i = 0; i < sendcounts[rank]; i++) {
        received_current_generation[i].chromosomes = malloc((c->size + 1) * sizeof(int));
    }

    /* auxiliary array of chromosomes used when receiving by all the workers */
    aux_chromosomes =  malloc((c->size + 1) * sizeof(int));
    

    /* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (int i = 0; i < generations_no; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        /* Scatering the current generation and sending the chromosomes */
        MPI_Scatterv(current_generation, sendcounts, displs, mpi_individual_type, received_current_generation, sendcounts[rank], mpi_individual_type, 0, MPI_COMM_WORLD);
        if (rank == ROOT) {
            int rank_count = 0;
            int local_pos = 0;
            for (int j = 0; j < population_size; j++) {
                if (local_pos >= sendcounts[rank_count]) {
                    local_pos = 0;
                    rank_count++;
                }
                MPI_Isend(current_generation[j].chromosomes, (c->size + 1), MPI_INT, rank_count, local_pos, MPI_COMM_WORLD, &request);
                local_pos++;
            }
        }
            
        for (int j = 0; j < sendcounts[rank]; j++) {
            MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            memcpy(received_current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
        }
        MPI_Barrier(MPI_COMM_WORLD);

        /* Calculating fitness in each worker */
        compute_generation_fitness_mpi_pthreads(received_current_generation, c, starting_point, sendcounts[rank], no_threads);
        MPI_Barrier(MPI_COMM_WORLD);

        /* Gathering in ROOT the individuals */
        MPI_Gatherv(received_current_generation, sendcounts[rank], mpi_individual_type, current_generation, sendcounts, displs, mpi_individual_type, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* send chromosomes from workers to ROOT */
        for (int j = 0; j < sendcounts[rank]; j++) {
            global_pos = global_pos_calc(rank, j, sendcounts, proc);
            MPI_Isend(received_current_generation[j].chromosomes, (c->size + 1), MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
        }

        /* reciving chromosomes from workers in ROOT */
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
            generate_random_numbers_mpi_pthreads(current_generation, c->size, population_size);

            /* adding the current generation and the next generation in the information for pthreads */
            for (int i = 0; i < no_threads; i++) {
                information[i].current_generation = current_generation;
                information[i].next_generation = next_generation;
            }

            /* creating the pthreads */
            int err;
            for (int i = 0; i < no_threads; i++) {
                err = pthread_create(&threads[i], NULL, (void *)pthreads_mutate, &information[i]);
                if (err) {
                    printf("Error when creating the thread\n");
                    exit(-1);
                }
            }

            /* joining the threads when mutating is completed */
            for (int i = 0; i < no_threads; i++) {
                int err = pthread_join(threads[i], NULL);

                if (err) {
                    printf("Eroare la asteptarea thread-ului %d\n", i);
                    exit(-1);
                }
            }

            /* Switch to new generation */
            auxiliary = current_generation;
            current_generation = next_generation;
            next_generation = auxiliary;
        }
        
        
    }
    

    if (rank == ROOT) {

        compute_generation_fitness_mpi_pthreads(current_generation, c, starting_point, population_size, no_threads);

        qsort(current_generation, population_size,
            sizeof(individual), compare_individuals_mpi);
        
        t2 = omp_get_wtime();
        printf("Total execution time = %lf\n", t2 - t1);
        print_result_individual_mpi(current_generation, c);
    }
    

	MPI_Finalize();
}
