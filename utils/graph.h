#pragma once

#include <stdio.h>
#include <stdlib.h>

/*
 * The structure that defines the graph
 * of cities and the costs of the edges
 */
typedef struct {
    int size;
    int **roads;
} cities;

void print_graph(cities *c) {
    int i, j;

    printf("%d\n", c->size);
    for (i = 0; i < c->size; i++) {
        printf("%d -> ", i + 1);
        for (j = 0; j < c->size; j++) {
            if (c->roads[i][j] != 0)
                printf("(%d, %d)", j, c->roads[i][j]);
        }
        printf("\n");
    }
}

cities* read_file(char *file) {
    FILE *f = fopen(file, "r");
    cities *cities_graph = malloc(sizeof(cities));
    int size, **c, i, j;

    fscanf(f, "%d", &size);
    c = malloc(size * sizeof(int*));
    for (i = 0; i < size; i++) {
        c[i] = malloc(size * sizeof(int));
        for (j = 0; j < size; j++) {
            fscanf(f, "%d", &c[i][j]);
        }
    }
    cities_graph->size = size;
    cities_graph->roads = c;
    return cities_graph;
}
