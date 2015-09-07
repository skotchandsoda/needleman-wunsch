// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// table.h - prototypes for table.c

#include <pthread.h>

#ifndef TABLE_H
#define TABLE_H

// cell_t, a type describing a cell in the table
typedef struct cell {
        int score;
        int diag;
        int left;
        int up;
        int match;
        int in_optimal_path;
        int processed;
        pthread_mutex_t score_mutex;
        pthread_cond_t processed_cv;
} cell_t;

// table_t, a type describing an MxN table of cells (cell_t)
typedef struct table {
        int M;
        int N;
        cell_t **cells;
        int greatest_abs_val;
} table_t;

// Allocate a MxN table of cells
table_t *alloc_table(int M, int N);

// Destroy a table of M cell_t pointers
void free_table(table_t *T, int multiple_threads);

// Initialize the score table
void init_table(table_t *T, int d, int multiple_threads);

// Print the score table
void print_table(table_t *T, char *s1, char *s2, int unicode);

#endif /* TABLE_H */
