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
 * Compare function used for sorting the individuals
 * from a certain generation
 */
int compare_individuals(const void *a, const void *b) {
    individual *individual_a = *(individual **) a;
    individual *individual_b = *(individual **) b;

    return individual_a->fitness - individual_b->fitness;
}


/*
 * Function that generates a series of random chromosomes
 * to serve as first generation
 */
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