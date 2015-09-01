// Copyright (C) 2015 Scott Cheloha.
// All rights reserved.

// cell.h -- prototypes for cell.c

#ifndef CELL_H
#define CELL_H

// A direction attribute for cells in the table
typedef enum {NONE, LEFT, UP, DIAG} direction;

// cell_t, a type describing a cell in the table
typedef struct cell {
    int score;
    direction ptr;
} cell_t;

// Print a cell in the table.  Optionally print the cell's direction
// next to its score.
void print_cell(cell_t *c, int print_direction);

#endif /* CELL_H */
