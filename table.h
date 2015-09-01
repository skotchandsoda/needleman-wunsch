// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// table.h -- prototypes for table.c

#include "cell.h"

#ifndef TABLE_H
#define TABLE_H

// Allocate a MxN table of cells
cell_t **alloc_table(int M, int N);

// Destroy a table of M cell_t pointers
void free_table(cell_t **table, int M);

// Initialize the score table
void init_table(cell_t **table, int M, int N, int d);

// Print the score table
void print_table(cell_t **table, int M, int N, char *s1, char *s2);

#endif /* TABLE_H */
