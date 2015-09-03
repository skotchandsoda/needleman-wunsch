// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// table.c -- routines for allocating, initializing, and printing
//            the scoring table in the Needleman Wunsch algorithm

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "format.h"
#include "table.h"

// A zeroed cell_t for when we first allocate the table
const cell_t DEFAULT_CELL = {0, 0, 0, 0, 0};

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
        // the cell_t arrays with the assignment to DEFAULT_CELL
        for (int i = 0; i < T->M; i++) {
                T->cells[i] = (cell_t *)malloc(T->N * sizeof(cell_t));
                if (T->cells[i] == NULL) {
                        perror("malloc failed");
                        exit(1);
                } else {
                        for (int j = 0; j < T->N; j++) {
                                T->cells[i][j] = DEFAULT_CELL;
                        }
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
                /* T->cells[i][0].up = 0; */
                /* T->cells[i][0].diag = 0; */
        }

        // Rest of leftmost column has score (-j)*d and UP direction
        for (int j = 1; j < T->N; j++) {
                T->cells[0][j].score = (-j) * d;
                T->cells[0][j].up = 1;
                /* T->cells[0][j].left = 0; */
                /* T->cells[0][j].diag = 0; */
        }
}

static void
print_top_row(table_t *T, char *s1, int score_col_width)
{
        printf("*    %*s", score_col_width, "*");
        for (int i = 0; i < T->M - 1; i++) {
                printf("    %*s%c", score_col_width-1, "", s1[i]);
        }

        printf("\n");
}

static void
print_row(table_t *T, int n, char *s1, char *s2, int score_col_width, int unicode)
{
        // Start with a space as a character placeholder
        printf(" ");

        // Print the row's header
        for (int i = 0; i < T->M; i++) {
                if (T->cells[i][n].diag == 1) {
                        if (T->cells[i][n].in_optimal_path == 1) {
                                if (s1[i] == s2[n]) {
                                        set_fmt(match_fmt);
                                } else {
                                        set_fmt(mismatch_fmt);
                                }
                        }
                        printf("  %s ", (unicode == 1 ? "\u2196" : "\\"));
                        reset_fmt();
                } else {
                        printf("    ");
                }

                if (T->cells[i][n].up == 1) {
                        // The unicode arrows are two bytes, so we add
                        // to the spacing offset when printing them
                    if (T->cells[i][n].in_optimal_path == 1) {
                        set_fmt(gap_fmt);
                    }
                        printf("%*s",
                               (unicode == 1 ? score_col_width + 2: score_col_width),
                               (unicode == 1 ? "\u2191" : "^"));
                        if (T->cells[i][n].in_optimal_path == 1) {
                            reset_fmt();
                        }
                } else {
                        printf("%*s", score_col_width, "");
                }
        }

        // Next, the row itself.  Start with either a '*' separator (if this
        // is the first row of numbers, i.e. n == 0) or a letter from the side
        // string (s2)
        printf("\n%c", (n == 0 ? '*' : s2[n-1]));

        // Now print the scores and left arrows
        for (int i = 0; i < T->M; i++) {
                if (T->cells[i][n].left == 1) {
                    if (T->cells[i][n].in_optimal_path == 1) {
                        set_fmt(gap_fmt);
                    }
                        printf("  %s ", (unicode == 1 ? "\u2190" : "<"));
                        if (T->cells[i][n].in_optimal_path == 1) {
                            reset_fmt();
                        }
                } else {
                        printf("    ");
                }
                if (T->cells[i][n].in_optimal_path == 1) {
                    set_fmt(bold_fmt);
                }
                printf("%+*d", score_col_width, T->cells[i][n].score);
                if (T->cells[i][n].in_optimal_path == 1) {
                    reset_fmt();
                }
        }

        printf("\n");
}

static int
width_needed_to_print(int x)
{
        int w = 0;
        do {
                x /= 10;
                w = w + 1;
        } while (x != 0);
        return w + 1; // add 1 to make room for a negative sign
}

void
print_table(table_t *T, char *s1, char *s2, int unicode)
{
        int score_col_width = width_needed_to_print(T->greatest_abs_val);
        print_top_row(T, s1, score_col_width);
        for (int i = 0; i < T->N; i++) {
                print_row(T, i, s1, s2, score_col_width, unicode);
        }
}
