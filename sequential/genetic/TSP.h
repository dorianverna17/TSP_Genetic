#include "../../utils/graph.h"

#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <limits.h>

typedef struct {
    int fitness;
    int *chromosomes;
    int position;
    int random_pos1;
    int random_pos2;
    int random_pos3;
    int random_pos4;
} individual;


/*
 * Function that computes the fitness of each individual
 * of a certain generation at some point in time.
 * Fitness is the cost of the journey between these cities
 */
void compute_generation_fitness(individual **generation, cities *c,
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

void generate_random_chromosomes(int *chromosomes, cities *c, int starting_point) {
    chromosomes[0] = starting_point;
    for (int i = 1; i < c->size; i++) {
        int city = rand() % c->size;

        while (check_chromosome(city, chromosomes, i))
            city = rand() % c->size;
        chromosomes[i] = city;
    }
    chromosomes[c->size] = starting_point;
}

/*
 * Function that generates two random numbers for
 * mutation swapping
 */
void generate_random_numbers(individual **gen, int size, int population_size) {
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
void minimum_chromosome_road(int *chromosomes, int start, cities *c) {
    int min, i, j, pos, city = start;

    chromosomes[0] = start;
    for (i = 0; i < c->size - 1; i++) {
        min = INT_MAX;
        for (j = 0; j < c->size; j++) {
            if (j != start && j != city &&
                c->roads[city][j] < min &&
                !check_chromosome(j, chromosomes, i)) {
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
void mutate_generation(individual **current_generation,
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
    generate_random_numbers(current_generation, c->size, population_size);
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
        generate_random_chromosomes(next_generation[i]->chromosomes, c, start);
    }
}

/*
 * Cities are taken as genes
 * Fitness score is the path length of all the cities mentioned
 */
void TSP_sequential_genetic(cities *c, int starting_point,
    int generations_no, int population_size) {
 
    /* Step 1: Creating initial population */

    individual **current_generation = malloc(population_size * sizeof(individual*));
    individual **next_generation = malloc(population_size * sizeof(individual*));
    individual **auxiliary;

    for (int i = 0; i < population_size; i++) {
        current_generation[i] = malloc(sizeof(individual));
        current_generation[i]->fitness = 0;
        current_generation[i]->position = i;
        current_generation[i]->chromosomes = malloc((c->size + 1) * sizeof(int));

        if (i != 0) {
            generate_random_chromosomes(current_generation[i]->chromosomes, c, starting_point);
        } else {
            minimum_chromosome_road(current_generation[i]->chromosomes, starting_point, c);
        }

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
        compute_generation_fitness(current_generation, c, starting_point, population_size);
        /* Step 3: Sort in order of fitnesses */
        qsort(current_generation, population_size,
            sizeof(individual*), compare_individuals);

        /* Step 4: Selecting the best genes and mutating */
        mutate_generation(current_generation, next_generation, c,
            starting_point, i, population_size);

        /* Step 5: Switch to new generation */
        auxiliary = current_generation;
        current_generation = next_generation;
        next_generation = auxiliary;
    }

    compute_generation_fitness(current_generation, c, starting_point, population_size);
    qsort(current_generation, population_size,
        sizeof(individual*), compare_individuals);

    print_result_individual(current_generation, c);
}
