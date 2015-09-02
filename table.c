// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// table.c -- routines for allocating, initializing, and printing
//            the scoring table in the Needleman Wunsch algorithm

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "table.h"

table_t *
alloc_table(int M, int N)
{
    // Allocate for the struct table
    table_t *T = (table_t *)malloc(sizeof(table_t));
    if (T == NULL) {
        perror("malloc failed");
        exit(1);
    }

    T->M = M;
    T->N = N;

    // Allocate top-level array for T->cells
    T->cells = (cell_t **)malloc(T->M * sizeof(cell_t *));
    if (T->cells == NULL) {
        perror("malloc failed");
        exit(1);
    }

    // Allocate subarrays for T->cells.  Note that we explicitly zero
    // the cell_t arrays with calloc(3).
    for (int i = 0; i < T->M; i++) {
      T->cells[i] = (cell_t *)calloc(T->N, sizeof(cell_t));
        if (T->cells[i] == NULL) {
            perror("malloc failed");
            exit(1);
        }
    }

    return T;
}

void
free_table(table_t *T)
{
    // Free each subarray of cells
    for (int i = 0; i < T->M; i++) {
        free(T->cells[i]);
    }

    // Free the top-level array of cell pointers
    free(T->cells);

    // Free the table_t itself
    free(T);
}

void
init_table(table_t *T, int d)
{
    // Initialize the table
    // 0,0 is 0 and has no direction
    T->cells[0][0].score = 0;


    // Rest of top row has score (-i)*d and LEFT direction
    for (int i = 1; i < T->M; i++) {
        T->cells[i][0].score = (-i) * d;
        T->cells[i][0].left = 1;
        T->cells[i][0].up = 0;
        T->cells[i][0].diag = 0;
    }

    // Rest of leftmost column has score (-j)*d and UP direction
    for (int j = 1; j < T->N; j++) {
        T->cells[0][j].score = (-j) * d;
        T->cells[0][j].up = 1;
        T->cells[0][j].left = 0;
        T->cells[0][j].diag = 0;
    }
}

static void
print_top_row(table_t *T, char *s1)
{
  printf("*     *");
  for (int i = 0; i < T->M - 1; i++) {
    printf("     %c", s1[i]);
  }

  printf("\n");
}

static void
print_row(table_t *T, int n, char *s2, int unicode)
{
  // Start with a '   ' separator
  printf(" ");

  // Print the row's header
  for (int i = 0; i < T->M; i++) {
    if (T->cells[i][n].diag == 1) {
      printf("  %s", (unicode == 1 ? "\u2B09" : "\\"));
    } else {
      printf("   ");
    }

    if (T->cells[i][n].up == 1) {
      printf("  %s", (unicode == 1 ? "\u2B06" : "\\"));
    } else {
      printf("   ");
    }
  }

  // Next, the row itself.  Start with either a '*' separator (if this
  // is the first row of numbers, i.e. n == 0) or a letter from the side
  // string (s2)
  printf("\n%c", (n == 0 ? '*' : s2[n-1]));

  // Now the scores and left separators
  for (int i = 0; i < T->M; i++) {
    if (T->cells[i][n].left == 1) {
      printf("  %s", (unicode == 1 ? "\u2B05" : "<"));
    } else {
      printf("   ");
    }
    printf("%3d", T->cells[i][n].score);
  }

  printf("\n");
}

void
print_table(table_t *T, char *s1, char *s2, int unicode)
{
  print_top_row(T, s1);
  for (int i = 0; i < T->N; i++) {
    print_row(T, i, s2, unicode);
  }
}
