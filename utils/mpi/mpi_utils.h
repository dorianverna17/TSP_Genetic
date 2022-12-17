#pragma once

#include <stdbool.h>
#include <string.h>
#include "../genetic_utils.h"
#include "mpi.h"
#include <stddef.h>

/*
 * Compare function used for sorting the individuals
 * from a certain generation
 */
int compare_individuals_mpi(const void *a, const void *b) {
    individual *individual_a = (individual *) a;
    individual *individual_b = (individual *) b;

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
void print_result_individual_mpi(individual *gen, cities *c) {
    int *journey = gen[0].chromosomes;

    printf("%d ", gen[0].fitness);
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

void print_generation_mpi(individual *gen, cities *c, int size) {
    for (int i = 0; i < size; i++) {
        print_individual_mpi(gen[i].chromosomes, c->size + 1, gen[i].fitness);
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
 * Defines custom MPI Datatype for individual so we can use them in MPI functions
*/
void define_MPI_individual_type(MPI_Datatype* mpi_individual) {
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

/*
 * Calculates the global position of an individual from a generation array
 * given the position in the chunk that a worker is using
*/
int global_pos_calc(int rank, int pos, int *sendcounts, int nrproc) {
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