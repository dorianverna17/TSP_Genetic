#include "../../utils/graph.h"

#include <stdbool.h>
#include <limits.h>

bool check_perm(int elem, int *perm, int size) {
    for (int i = 0; i < size; i++) {
        if (elem == perm[i]) {
            return false;
        }
    }
    return true;
}

void print_res(int *perm, int step, int start) {
    for (int i = 0; i < step; i++) {
        printf("%d -> ", perm[i] + 1);
    }
    printf("%d\n", start + 1);
}

void tsp_backtrack(int *perm, cities *c, int step, int *cost, int *res, int aux_cost, int start) {
    if (step == c->size) {
        if (c->roads[perm[step - 1]][start] == 0) {
            return;
        } else {
            aux_cost += c->roads[perm[step - 1]][start];
        }
        if (*cost > aux_cost) {
            /* copy permutation into the result array */
            for (int i = 0; i < c->size; i++) {
                res[i] = perm[i];
            }
            *cost = aux_cost;
        }
        return;
    }
    for (int i = 0; i < c->size; i++) {
        if (check_perm(i, perm, step)) {
            perm[step] = i;
        } else {
            continue;
        }
        if (c->roads[perm[step - 1]][perm[step]] == 0)
            continue;
        if (step > 0) {
            tsp_backtrack(perm, c, step + 1, cost, res,
                aux_cost + c->roads[perm[step - 1]][perm[step]], start);
        } else {
            tsp_backtrack(perm, c, step + 1, cost, res, aux_cost, start);
        }
        perm[step] = -1;
    }
}

void TSP_sequential_naive(cities *c, int starting_point) {
    /* step 1: generate all (n - 1)! permutaions of cities */
    int *perm = malloc(c->size * sizeof(int));
    int *result = malloc(c->size * sizeof(int));
    int cost = INT_MAX;

    perm[0] = starting_point;
    for (int i = 1; i < c->size; i++) {
        perm[i] = -1;
    }

    tsp_backtrack(perm, c, 1, &cost, result, 0, starting_point);

    print_res(result, c->size, starting_point);
    printf("%d\n", cost);
}