#include "../../utils/graph.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>

typedef struct {
    int fitness;
    int *chromosomes;
    int position;
} individual;

/*
 * Rotate the position of the cities using the
 * index given as second parameter
 */
void generate_custom_chromosomes(int *chromosomes, int position,
    int start, int size) {

    if (position == size - 2) {
        for (int i = size - 2; i >= 0; i--) {
            chromosomes[size - 1 - i] = i;
        }
        chromosomes[0] = start;
        return;
    }

    chromosomes[0] = start;

    for (int i = 1; i <= position; i++) {
        chromosomes[i] = size - 2 - position + i;
    }
    for (int i = position + 1; i < size - 1; i++) {
        chromosomes[i] = i - position;
    }

    chromosomes[size - 1] = start;
}

/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness(individual **generation, cities *c, int start) {
    int cost;
    int *chromosomes;

    for (int i = 0; i < c->size; i++) {
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
int compare_individuals(const void *a, const void *b) {
    individual *individual_a = *(individual **) a;
    individual *individual_b = *(individual **) b;

    return individual_a->fitness - individual_b->fitness;
}

/*
 * Function that checks if a certain chromosome is
 * already in the array of chromosomes
 */
bool check_chromosome(int chromosome, int *chromosomes, int size) {
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
void print_result_individual(individual **gen, cities *c) {
    int *journey = gen[0]->chromosomes;

    printf("%d ", gen[0]->fitness);
    for (int i = 0; i < c->size; i++) {
        printf("%d -> ", journey[i] + 1);
    }
    printf("%d\n", journey[c->size] + 1);
}

void print_individual(int *chromosomes, int size, int fitness) {
    printf("%d ", fitness);

    for (int i = 0; i < size - 1; i++) {
        printf("%d -> ", chromosomes[i]);
    }
    printf("%d\n", chromosomes[size - 1]);
}

void print_generation(individual **gen, cities *c) {
    for (int i = 0; i < c->size; i++) {
        print_individual(gen[i]->chromosomes, c->size + 1, gen[i]->fitness);
    }
    printf("\n");
}

void generate_random_chromosomes(int *chromosomes, cities *c) {
    for (int i = 1; i < c->size; i++) {
        int city = rand() % c->size;

        while (check_chromosome(city, chromosomes, i))
            city = rand() % c->size;
        chromosomes[i] = city;
    }
}

/*
 * Function that populates the next generation with the
 * first 40% individuals of the current generation, mutates
 * the first 30% and then crossover the first 30%
 */
void mutate_generation(individual **current_generation,
    individual **next_generation, cities *c, int start, int gen_no) {
    int i, j, aux;
    int current_index;

    /* keep the first 40% */
    int count_best = (c->size * 4) / 10;
    for (i = 0; i < count_best; i++) {
        memcpy(next_generation[i]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[i]->fitness = 0;
    }


    current_index = count_best;

    /* 
     * mutate the first 30% - swap the position of each
     * two consecutive chromosomes of an individual
     */
    int count = (c->size * 2) / 10;
    for (i = 0; i < count; i++) {
        memcpy(next_generation[current_index]->chromosomes,
            current_generation[current_index]->chromosomes, (c->size + 1) * sizeof(int));
        next_generation[current_index]->fitness = 1;
        for (j = 1; j < c->size - 1; j+=2) {
            aux = next_generation[current_index]->chromosomes[j];
            next_generation[current_index]->chromosomes[j] =
                next_generation[current_index]->chromosomes[j + 1];
            next_generation[current_index]->chromosomes[j + 1] = aux;
        }
        current_index += 1;
    }

    /* the rest of them are being mutated this way */
    for (i = 0; i < count; i++) {
        memcpy(next_generation[current_index]->chromosomes,
            current_generation[i]->chromosomes, (c->size + 1) * sizeof(int));
        for (int j = 1; j < c->size - 2; j+=3) {
            aux = next_generation[current_index]->chromosomes[j];
            next_generation[current_index]->chromosomes[j] =
                next_generation[current_index]->chromosomes[j + 2];
            next_generation[current_index]->chromosomes[j + 2] = aux;
        }
        next_generation[current_index]->fitness = 3;
        current_index += 1;
    }

    /* the remaining ones are being generated randomly */
    for (i = current_index; i < c->size; i++) {
        generate_random_chromosomes(next_generation[i]->chromosomes, c);
        next_generation[i]->fitness = 4;
    }

}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 * 
 */
void TSP_sequential_genetic(cities *c, int starting_point, int generations_no) {
    srand(time(NULL));
 
    /* Step 1: Creating initial population */

    individual **current_generation = malloc(c->size * sizeof(individual*));
    individual **next_generation = malloc(c->size * sizeof(individual*));
    individual **auxiliary;

    for (int i = 0; i < c->size; i++) {
        current_generation[i] = malloc(sizeof(individual));
        current_generation[i]->fitness = 0;
        current_generation[i]->position = i;
        current_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        generate_custom_chromosomes(current_generation[i]->chromosomes,
            i, starting_point, c->size + 1);

        next_generation[i] = malloc(sizeof(individual));
        next_generation[i]->fitness = 0;
        next_generation[i]->position = i;
        next_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));
    
        for (int j = 0; j < c->size + 1; j++) {
            next_generation[i]->chromosomes[j] = 0;
        }
    }

    /* 
     * Complete the other steps until we reach the desired
     * number of generations
     */
    for (int i = 0; i < generations_no; i++) {
        /* Step 2: Calculating fitness */
        compute_generation_fitness(current_generation, c, starting_point);
        /* Step 3: Sort in order of fitnesses */
        qsort(current_generation, c->size,
            sizeof(individual*), compare_individuals);

        // print_generation(current_generation, c);

        /* Step 4: Selecting the best genes and mutating */
        mutate_generation(current_generation, next_generation, c, starting_point, i);
        
        // print_generation(next_generation, c);

        /* Step 5: Switch to new generation */
        auxiliary = current_generation;
        current_generation = next_generation;
        next_generation = auxiliary;
    }

    compute_generation_fitness(current_generation, c, starting_point);
    qsort(current_generation, c->size,
        sizeof(individual*), compare_individuals);

    print_result_individual(current_generation, c);
}