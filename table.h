// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// table.h -- prototypes for table.c

#ifndef TABLE_H
#define TABLE_H

// cell_t, a type describing a cell in the table
typedef struct cell {
        int score;
        int diag;
        int left;
        int up;
        int in_optimal_path;
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
void free_table(table_t *T);

// Initialize the score table
void init_table(table_t *T, int d);

// Print the score table
void print_table(table_t *T, char *s1, char *s2, int unicode);

#endif /* TABLE_H */
