// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// table.c - routines for allocating, initializing, and printing
//           the scoring table in the Needleman Wunsch algorithm

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "format.h"
#include "table.h"

#define UNICODE_LEFTWARDS_ARROW "\u2190"
#define UNICODE_UPWARDS_ARROW "\u2191"
#define UNICODE_NORTH_WEST_ARROW "\u2196"

// A zeroed cell_t for when we first allocate the table
const cell_t DEFAULT_CELL = {0, 0, 0, 0, 0, 0, 0};

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
        T->cells[0][0].processed = 1;


        // Rest of top row has score (-i)*d and LEFT direction
        for (int i = 1; i < T->M; i++) {
                T->cells[i][0].score = (-i) * d;
                T->cells[i][0].left = 1;
                T->cells[i][0].processed = 1;
        }

        // Rest of leftmost column has score (-j)*d and UP direction
        for (int j = 1; j < T->N; j++) {
                T->cells[0][j].score = (-j) * d;
                T->cells[0][j].up = 1;
                T->cells[0][j].processed = 1;
        }
}

typedef enum {left, up, diag} arrow_t;

static void
print_arrow(arrow_t a, int optimal_path, int col_width, int col, int row, char *s1, char *s2, int unicode)
{
    switch (a) {
    case left:
            if (optimal_path == 1) {
                    set_fmt(gap_arrow_fmt);
            }
            printf("  %s ", (unicode == 1 ? UNICODE_LEFTWARDS_ARROW : "<"));
            break;
    case up:
            if (optimal_path == 1) {
                    set_fmt(gap_arrow_fmt);
            }
            printf("%*s",
                   (unicode == 1 ? col_width + 2: col_width),
                   (unicode == 1 ? UNICODE_UPWARDS_ARROW : "^"));
            break;
    case diag:
            if (optimal_path == 1) {
                    // If the characters match, use the matching arrow
                    // format.  If the characters do not match, use the
                    // mismatching arrow format.
                    fmt_t f = s1[col-1] == s2[row-1] ? match_arrow_fmt : mismatch_arrow_fmt;
                    set_fmt(f);
            }
            printf("  %s ", (unicode == 1 ? UNICODE_NORTH_WEST_ARROW : "\\"));
            break;
    default:
            fprintf(stderr, "the impossible has happened; giving up\n");
            exit(1);
            break;
    }

    if (optimal_path == 1) {
            reset_fmt();
    }
}

static void
print_directional_row(table_t *T, int row, char *s1, char *s2, int col_width, int unicode)
{
        // Start with a space as a character placeholder
        printf(" ");

        // Print the row's directional arrows
        for (int col = 0; col < T->M; col++) {
                int optimal_path = T->cells[col][row].in_optimal_path;

                // Print diagonal arrow if applicable
                if (T->cells[col][row].diag == 1) {
                        print_arrow(diag, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("    ");
                }

                // Print up arrow if applicable
                if (T->cells[col][row].up == 1) {
                        print_arrow(up, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("%*s", col_width, "");
                }
        }
        printf("\n");
}

static void
print_score_row(table_t *T, int row, char *s1, char *s2, int col_width, int unicode)
{
        // Start with either a '*' separator (if this
        // is the first row of numbers, i.e. n == 0) or a letter from the side
        // string (s2)
        set_fmt(side_string_fmt);
        printf("%c", (row == 0 ? '*' : s2[row-1]));
        reset_fmt();

        // Now print the scores and left arrows
        for (int col = 0; col < T->M; col++) {
                int optimal_path = T->cells[col][row].in_optimal_path;

                // Print left arrow if applicable
                if (T->cells[col][row].left == 1) {
                        print_arrow(left, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("    ");
                }

                // Print the cell's score
                if (optimal_path == 1) {
                        set_fmt(opt_path_fmt);
                }
                printf("%+*d", col_width, T->cells[col][row].score);
                if (optimal_path == 1) {
                        reset_fmt();
                }
        }
        printf("\n");
}

static void
print_table_row(table_t *T, int row, char *s1, char *s2, int col_width, int unicode)
{
        print_directional_row(T, row, s1, s2, col_width, unicode);
        print_score_row(T, row, s1, s2, col_width, unicode);
}

static void
print_top_string(table_t *T, char *s1, int col_width)
{
        set_fmt(top_string_fmt);
        printf("*    %*s", col_width, "*");
        for (int i = 0; i < T->M - 1; i++) {
                printf("    %*s%c", col_width-1, "", s1[i]);
        }

        printf("\n");
}

static int
width_needed_to_print_integer(int x)
{
        int w = 0;
        do {
                x /= 10;
                w = w + 1;
        } while (x != 0);
        return w + 1; // add 1 to make room for a positive/negative sign
}

void
print_table(table_t *T, char *s1, char *s2, int unicode)
{
        int col_width = width_needed_to_print_integer(T->greatest_abs_val);

        // Print the "top string"
        print_top_string(T, s1, col_width);

        // Print the rest of the rows, bordered on the left by the "side string"
        for (int i = 0; i < T->N; i++) {
                print_table_row(T, i, s1, s2, col_width, unicode);
        }
}
