#include "../../utils/graph.h"
#include "../../utils/genetic_utils.h"
#include "mpi.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <math.h>
#include <omp.h>

#define ROOT 0

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness_mpi_omp(individual *generation, cities *c,
    int start, int population_size, int no_threads) {
    int cost, start_index, end_index, i, thread_id, j;
    int *chromosomes;

    omp_set_num_threads(no_threads);

    #pragma omp parallel private(thread_id, start_index, end_index, i, cost, j, chromosomes) shared(generation, c, no_threads, population_size)
	{
		thread_id = omp_get_thread_num();
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
        #pragma omp barrier
    }
}

/*
 * Compare function used for sorting the individuals
 * from a certain generation
 */
int compare_individuals_mpi_omp(const void *a, const void *b) {
    individual *individual_a = (individual *) a;
    individual *individual_b = (individual *) b;

    return individual_a->fitness - individual_b->fitness;
}

/*
 * Function that checks if a certain chromosome is
 * already in the array of chromosomes
 */
bool check_chromosome_mpi_omp(int chromosome, int *chromosomes, int size) {
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
void print_result_individual_mpi_omp(individual *gen, cities *c) {
    int *journey = gen[0].chromosomes;

    printf("%d ", gen[0].fitness);
    for (int i = 0; i < c->size; i++) {
        printf("%d -> ", journey[i] + 1);
    }
    printf("%d\n", journey[c->size] + 1);
}

void print_individual_mpi_omp(int *chromosomes, int size, int fitness) {
    printf("%d ", fitness);

    for (int i = 0; i < size - 1; i++) {
        printf("%d -> ", chromosomes[i]);
    }
    printf("%d\n", chromosomes[size - 1]);
}

void print_generation_mpi_omp(individual *gen, cities *c, int size) {
    for (int i = 0; i < size; i++) {
        print_individual_mpi_omp(gen[i].chromosomes, c->size + 1, gen[i].fitness);
    }
    printf("\n");
}

void generate_random_chromosomes_mpi_omp(int *chromosomes, cities *c, int starting_point) {
    chromosomes[0] = starting_point;
    for (int i = 1; i < c->size; i++) {
        int city = rand() % c->size;

        while (check_chromosome_mpi_omp(city, chromosomes, i))
            city = rand() % c->size;
        chromosomes[i] = city;
    }
    chromosomes[c->size] = starting_point;
}

/*
 * Function that generates two random numbers for
 * mutation swapping
 */
void generate_random_numbers_mpi_omp(individual *gen, int size, int population_size) {

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

/*
 * Function that constructs an individual by choosing the
 * minimum distance from each city when that particular
 * city is visited by the travelling salesman
 */
void minimum_chromosome_road_mpi_omp(int *chromosomes, int start, cities *c) {
    int min, i, j, pos, city = start;

    chromosomes[0] = start;
    for (i = 0; i < c->size - 1; i++) {
        min = INT_MAX;
        for (j = 0; j < c->size; j++) {
            if (j != start && j != city &&
                c->roads[city][j] < min &&
                !check_chromosome_mpi_omp(j, chromosomes, i)) {
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
void mutate_generation_mpi_omp(individual *current_generation,
    individual *next_generation, cities *c, int start, int gen_no, int population_size, int no_threads) {
    int i, j, aux, thread_id, current_index, start_index, end_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;

    omp_set_num_threads(no_threads);

    #pragma omp parallel private(thread_id, start_index, end_index, i, current_index, aux) shared(current_generation, next_generation, c, no_threads, population_size, count_best)
	{
        thread_id = omp_get_thread_num();
        start_index = thread_id * (double) population_size / no_threads;
        end_index = min((thread_id + 1) * (double) population_size / no_threads, population_size);

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

        current_index = i;

        // let's mutate the rest of them
        for (i = max(start_index, current_index); i < min(current_index + count_best, end_index); i++) {
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
}

/*
 * Defines custom MPI Datatype for individual so we can use them in MPI functions
*/
void define_MPI_individual_type_mpi_omp(MPI_Datatype* mpi_individual) {
    MPI_Datatype oldtypes[6];
    int blockcounts[6];
    MPI_Aint offsets[6];
    MPI_Status status;

    offsets[0] = offsetof(individual, fitness);
    oldtypes[0] = MPI_INT;
    blockcounts[0] = 1;

    offsets[1] = offsetof(individual, position);
    oldtypes[1] = MPI_INT;
    blockcounts[1] = 1;

    offsets[2] = offsetof(individual, random_pos1);
    oldtypes[2] = MPI_INT;
    blockcounts[2] = 1;

    offsets[3] = offsetof(individual, random_pos2);
    oldtypes[3] = MPI_INT;
    blockcounts[3] = 1;

    offsets[4] = offsetof(individual, random_pos3);
    oldtypes[4] = MPI_INT;
    blockcounts[4] = 1;

    offsets[5] = offsetof(individual, random_pos4);
    oldtypes[5] = MPI_INT;
    blockcounts[5] = 1;

    MPI_Datatype mpi_struct;
    MPI_Type_create_struct(6, blockcounts, offsets, oldtypes, &mpi_struct);
    MPI_Type_create_resized(mpi_struct, 0, sizeof(individual), mpi_individual);
    MPI_Type_commit(mpi_individual);
}

int global_pos_calc_mpi_omp(int rank, int pos, int *sendcounts, int nrproc) {
    int global_pos = 0;
    for (int i = 0; i < nrproc; i++) {
        if (rank == i) {
            global_pos += pos;
            break;
        } else {
            global_pos += sendcounts[i];
        }
    }
    return global_pos;
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_parallel_mpi_omp(cities *c, int starting_point, int generations_no, int population_size, int no_threads) {

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
    define_MPI_individual_type_mpi_omp(&mpi_individual_type);

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
                minimum_chromosome_road_mpi_omp(current_generation[i].chromosomes, starting_point, c);
            } else {
                generate_random_chromosomes_mpi_omp(current_generation[i].chromosomes, c, starting_point);
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
        compute_generation_fitness_mpi_omp(received_current_generation, c, starting_point, sendcounts[rank], no_threads);
        MPI_Barrier(MPI_COMM_WORLD);

        /* Gathering in ROOT the individuals */
        MPI_Gatherv(received_current_generation, sendcounts[rank], mpi_individual_type, current_generation, sendcounts, displs, mpi_individual_type, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        /* send chromosomes from workers to ROOT */
        for (int j = 0; j < sendcounts[rank]; j++) {
            global_pos = global_pos_calc_mpi_omp(rank, j, sendcounts, proc);
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
                sizeof(individual), compare_individuals_mpi_omp); 
            generate_random_numbers_mpi_omp(current_generation, c->size, population_size);

            /* Step 4: Selecting the best genes and mutating */
            mutate_generation_mpi_omp(current_generation, next_generation, c,
                starting_point, i, population_size, no_threads);

            /* Step 5: Switch to new generation */
            auxiliary = current_generation;
            current_generation = next_generation;
            next_generation = auxiliary;
        }
        
        
    }
    

    if (rank == ROOT) {

        compute_generation_fitness_mpi_omp(current_generation, c, starting_point, population_size, no_threads);

        qsort(current_generation, population_size,
            sizeof(individual), compare_individuals_mpi_omp);
        
        t2 = omp_get_wtime();
        printf("Total execution time = %lf\n", t2 - t1);
        print_result_individual_mpi_omp(current_generation, c);
    }
    

	MPI_Finalize();
}
