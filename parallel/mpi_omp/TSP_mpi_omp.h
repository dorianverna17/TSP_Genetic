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
    int start, int population_size) {
    int cost;
    int *chromosomes;


    for (int i = 0; i < population_size; i++) {
        cost = 0;
        chromosomes = generation[i].chromosomes;

        for (int j = 0; j < c->size; j++) {
            cost += c->roads[chromosomes[j]][chromosomes[j + 1]];
        }
        generation[i].fitness = cost;
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
        printf("i = %d ", i);
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
    individual *next_generation, cities *c, int start, int gen_no, int population_size) {
    int i, j, aux;
    int current_index;

    /* keep the first 20% */
    int count_best = (population_size * 3) / 10;
    for (i = 0; i < count_best; i++) {
        memcpy(next_generation[i].chromosomes,
            current_generation[i].chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i].fitness = 0;
    }

    current_index = count_best;

    // let's mutate the rest of them
    generate_random_numbers_mpi_omp(current_generation, c->size, population_size);
    for (i = current_index; i < current_index + count_best; i++) {
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
    for (i = current_index; i < current_index + count_best; i++) {
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

    for (i = current_index; i < population_size; i++) {
        generate_random_chromosomes_mpi_omp(next_generation[i].chromosomes, c, start);
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
void TSP_parallel_mpi_omp(cities *c, int starting_point, int generations_no, int population_size) {

	int rank, proc, global_pos;
    double t1, t2; 
    
    individual *current_generation = NULL;
	individual *next_generation = NULL;
	individual *auxiliary = NULL;
    individual *received_current_generation = NULL;
    individual *received_next_generation = NULL;
    int *aux_chromosomes;

    int provided;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
    MPI_Request request = MPI_REQUEST_NULL;
    MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    t1 = MPI_Wtime();

    /* Step 1: Creating initial population */
    if (rank == ROOT) {
        current_generation = malloc(population_size * sizeof(individual));
        next_generation = malloc(population_size * sizeof(individual));
    }
        // received_current_generation = malloc(population_size * sizeof(individual));

	/* Create MPI Data_type for individuals*/
    MPI_Datatype mpi_individual_type;
    define_MPI_individual_type_mpi_omp(&mpi_individual_type);

	if (rank == ROOT) {
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
    
    /* the arrays all workers will use*/
    received_current_generation = (individual*)malloc(sendcounts[rank] * sizeof(individual));
    received_next_generation = (individual*)malloc(sendcounts[rank]*sizeof(individual));
    aux_chromosomes =  malloc((c->size + 1) * sizeof(int));

    MPI_Scatterv(current_generation, sendcounts, displs, mpi_individual_type, received_current_generation, sendcounts[rank], mpi_individual_type, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


    /*
     * Allocate memory for the chromosomes of each generation and initialize them
    */
    for (int i = 0; i < sendcounts[rank]; i++) {
        received_current_generation[i].chromosomes = malloc((c->size + 1) * sizeof(int));
        received_next_generation[i].chromosomes = calloc((c->size + 1), sizeof(int));
        if (rank == 0 && i == 0) {
            minimum_chromosome_road_mpi_omp(received_current_generation[i].chromosomes, starting_point, c);
        } else {
            generate_random_chromosomes_mpi_omp(received_current_generation[i].chromosomes, c, starting_point);
        }
    }
    

    /* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (int i = 0; i < generations_no; i++) {
        printf("----RANKUL %d la generatia %d-----\n", rank, i);
        /* Step 2: Calculating fitness */
        MPI_Barrier(MPI_COMM_WORLD);
        global_pos = global_pos_calc_mpi_omp(rank, 0, sendcounts, proc);
        compute_generation_fitness_mpi_omp(received_current_generation, c, starting_point, sendcounts[rank]);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(received_current_generation, sendcounts[rank], mpi_individual_type, current_generation, sendcounts, displs, mpi_individual_type, 0, MPI_COMM_WORLD);

        /* send chromosomes from workers to ROOT */
        for (int j = 0; j < sendcounts[rank]; j++) {
            global_pos = global_pos_calc_mpi_omp(rank, j, sendcounts, proc);
            MPI_Isend(received_current_generation[j].chromosomes, (c->size + 1), MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
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
            /* Step 3: Sort in order of fitnesses */
            qsort(current_generation, population_size,
                sizeof(individual), compare_individuals_mpi_omp); 
            
            // mutate cu omp -> scapam de bcast la random 
            // scatter - compute (cu omp in compute) - gather + isend cromozomi - qsort secv - mutate cu omp de pe un singur thread

            global_pos = global_pos_calc_mpi_omp(rank, 0, sendcounts, proc);

            /* Step 4: Selecting the best genes and mutating */
            mutate_generation_mpi_omp(current_generation, next_generation, c,
                starting_point, i, population_size);
        
            /* Step 5: Switch to new generation */
            auxiliary = current_generation;
            current_generation = next_generation;
            next_generation = auxiliary;
        }

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
                MPI_Isend(current_generation[j].chromosomes, (c->size + 1), MPI_INT, ROOT, local_pos, MPI_COMM_WORLD, &request);
                local_pos++;
            }
        }
            
        for (int j = 0; j < sendcounts[rank]; j++) {
            MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            memcpy(received_current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
        }
    }
    

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(received_current_generation, sendcounts[rank], mpi_individual_type, current_generation, sendcounts, displs, mpi_individual_type, 0, MPI_COMM_WORLD);

    /* send chromosomes from workers to ROOT */
    for (int j = 0; j < sendcounts[rank]; j++) {
        global_pos = global_pos_calc_mpi_omp(rank, j, sendcounts, proc);
        MPI_Isend(received_current_generation[global_pos].chromosomes, (c->size + 1), MPI_INT, ROOT, global_pos, MPI_COMM_WORLD, &request);
    
    }

    /* receiving the chromosomes */
    if (rank == ROOT) {
        for (int j = 0; j < population_size; j++) {
            MPI_Irecv(aux_chromosomes, c->size + 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            MPI_Wait(&request, &status);
            memcpy(received_current_generation[status.MPI_TAG].chromosomes, aux_chromosomes, (c->size + 1) * sizeof(int));
        }
    }

    if (rank == ROOT) {

        compute_generation_fitness_mpi_omp(current_generation, c, starting_point, population_size);

        qsort(current_generation, population_size,
            sizeof(individual), compare_individuals_mpi_omp);
        
        t2 = MPI_Wtime();
        printf("Total execution time = %lf\n", t2 - t1);
        print_result_individual_mpi_omp(current_generation, c);
    }
    

	MPI_Finalize();
}
