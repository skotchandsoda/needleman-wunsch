// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// table.c -- routines for allocating, initializing, and printing
//            the scoring table in the Needleman Wunsch algorithm

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "cell.h"

cell_t **
alloc_table(int M, int N)
{
    // Allocate for table top-level array
    cell_t **table = (cell_t **)malloc(M * sizeof(cell_t *));
    if (table == NULL) {
        perror("malloc failed");
        exit(1);
    }

    // Allocate subarrays
    for (int i = 0; i < M; i++) {
        table[i] = (cell_t *)malloc(N * sizeof(cell_t));
        if (table[i] == NULL) {
            perror("malloc failed");
            exit(1);
        }
    }

    return table;
}

void
free_table(cell_t **table, int M)
{
    // Free each subarray of cells
    for (int i = 0; i < M; i++) {
        free(table[i]);
    }

    // Free the top-level array of cell pointers
    free(table);
}

void
init_table(cell_t **table, int M, int N, int d)
{
    // Initialize the table
    // 0,0 is 0 and has no direction
    table[0][0].score = 0;
    table[0][0].ptr = NONE;

    // Rest of top row has score (-i)*d and LEFT direction
    for (int i = 1; i < M; i++) {
        table[i][0].score = (-i) * d;
        table[i][0].ptr = LEFT;
    }

    // Rest of leftmost column has score (-j)*d and UP direction
    for (int j = 1; j < N; j++) {
        table[0][j].score = (-j) * d;
        table[0][j].ptr = UP;
    }
}

void
print_table(cell_t **table, int M, int N, char *s1, char *s2)
{
    // print s1 padded with two '*' placeholders
    printf("* *");
    for (int i = 0; i < M; i++) {
        printf(" %c", s1[i]);
    }

    // Print first row of numbers preceded by '*' placeholder
    printf("\n*");
    for (int i = 0; i < M; i++) {
        printf(" ");
        print_cell(&table[i][0], 1);
    }

    // Print the rest of the rows, all preceded by a letter from s2
    for (int j = 1; j < N; j++) {
        printf("\n%c", s2[j-1]);
        for (int i = 0; i < M; i++) {
            printf(" ");
            print_cell(&table[i][j], 1);
        }
    }
    printf("\n");
}
