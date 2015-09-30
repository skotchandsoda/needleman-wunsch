/*-
 * Copyright (c) 2015, Scott Cheloha.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the copyright holder nor the names of its
 *      contributors may be used to endorse or promote products derived from
 *      this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * score-table.c - Routines for allocating, initializing, and printing
 *                 the scoring table in the Needleman Wunsch algorithm.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "dbg.h"
#include "format.h"
#include "score-table.h"
#include "walk-table.h"

#define UNICODE_LEFTWARDS_ARROW  "\u2190"
#define UNICODE_UPWARDS_ARROW    "\u2191"
#define UNICODE_NORTH_WEST_ARROW "\u2196"

score_table_t *
alloc_score_table(int M, int N)
{
        /* Allocate for the scores table */
        score_table_t *S = (score_table_t *)malloc(sizeof(score_table_t));
        check(NULL != S, "malloc failed");

        S->M = M;
        S->N = N;

        /* Allocate top-level array for S->cells */
        S->cells = (score_table_cell_t **)malloc(M * sizeof(score_table_cell_t *));
        check(NULL != S->cells, "malloc failed");

        /* Allocate subarrays for S->cells */
        for (int i = 0; i < M; i++) {
                S->cells[i] = (score_table_cell_t *)calloc(N, sizeof(score_table_cell_t));
                check(NULL != S->cells[i], "malloc failed");
        }

        return S;
}

void
free_score_table(score_table_t *S, unsigned int nthreads)
{
        /* Destroy mutex and conditional variable objects */
        int res;
        if (nthreads > 1) {
                for (int i = 0; i < S->M; i++) {
                        for (int j = 0; j < S->N; j++) {
                                res = pthread_mutex_destroy(&S->cells[i][j].score_mutex);
                                check(0 == res, "pthread_mutex_destroy failed");
                                res = pthread_cond_destroy(&S->cells[i][j].processed_cv);
                                check(0 == res, "pthread_cond_destroy failed");
                        }
                }
        }

        /* Free each subarray (column) of cells */
        for (int i = 0; i < S->M; i++) {
                free(S->cells[i]);
        }

        /* Free the top-level array of cell pointers */
        free(S->cells);

        /* Free the scores table itself */
        free(S);
}

/* void */
/* init_score_table(score_table_t *S, int d, int multiple_threads) */
/* { */
/*         /\* Initialize all mutex and condition variable objects *\/ */
/*         int res; */
/*         if (1 == multiple_threads) { */
/*                 for (int i = 0; i < S->M; i++) { */
/*                         for (int j = 0; j < S->N; j++) { */
/*                                 res = pthread_mutex_init(&S->cells[i][j].score_mutex, NULL); */
/*                                 check(0 == res, "pthread_mutex_init failed"); */
/*                                 res = pthread_cond_init (&S->cells[i][j].processed_cv, NULL); */
/*                                 check(0 == res, "pthread_cond_init failed"); */
/*                         } */
/*                 } */
/*         } */

/*         /\* Initialize the largest value *\/ */
/*         S->greatest_abs_val = 0; */

/*         /\* Initialize the table.  0,0 has score 0 and has no optimal */
/*            direction *\/ */
/*         S->cells[0][0].score = 0; */
/*         S->cells[0][0].processed = 1; */
/*         S->cells[0][0].up_done = 1; */
/*         S->cells[0][0].left_done = 1; */
/*         S->cells[0][0].diag_done = 1; */

/*         // Rest of top row has score i * (-d) and LEFT direction */
/*         for (int i = 1; i < S->M; i++) { */
/*                 S->cells[i][0].score = i * (-d); */
/*                 S->cells[i][0].left = 1; */
/*                 S->cells[i][0].up_done = 1; */
/*                 S->cells[i][0].diag_done = 1; */
/*                 S->cells[i][0].processed = 1; */
/*         } */

/*         /\* Rest of leftmost (side) column has score j * (-d) and UP */
/*          * direction *\/ */
/*         for (int j = 1; j < S->N; j++) { */
/*                 S->cells[0][j].score = j * (-d); */
/*                 S->cells[0][j].up = 1; */
/*                 S->cells[0][j].left_done = 1; */
/*                 S->cells[0][j].diag_done = 1; */
/*                 S->cells[0][j].processed = 1; */
/*         } */

/*         S->branch_count = 0; */
/*         res = pthread_rwlock_init(&(S->branch_count_rwlock), NULL); */
/*         check(0 == res, "pthread_rwlock_init failed"); */
/* } */


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
            unreachable();
            break;
    }

    if (optimal_path == 1) {
            reset_fmt();
    }
}

static void
print_directional_row(walk_table_t *W,
                      int row,
                      char *s1,
                      char *s2,
                      int col_width,
                      int unicode)
{
        // Start with a space as a character placeholder
        printf(" ");

        // Print the row's directional arrows
        for (int col = 0; col < W->M; col++) {
                int optimal_path = W->cells[col][row].in_optimal_path;

                /* Print diagonal arrow if applicable */
                if (W->cells[col][row].diag == 1) {
                        print_arrow(diag, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("    ");
                }

                /* Print up arrow if applicable */
                if (W->cells[col][row].up == 1) {
                        print_arrow(up, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("%*s", col_width, "");
                }
        }
        printf("\n");
}

static void
print_score_row(score_table_t *S,
                walk_table_t *W,
                int row,
                char *s1,
                char *s2,
                int col_width,
                int unicode)
{
        /* Start with either a '-' separator (if this
           is the first row of numbers, i.e. n == 0) or a letter from the side
           string (s2) */
        set_fmt(side_string_fmt);
        printf("%c", (row == 0 ? '-' : s2[row-1]));
        reset_fmt();

        // Now print the scores and left arrows
        for (int col = 0; col < S->M; col++) {
                int optimal_path = W->cells[col][row].in_optimal_path;

                /* Print left arrow if applicable */
                if (W->cells[col][row].left == 1) {
                        print_arrow(left, optimal_path, col_width, col, row, s1, s2, unicode);
                } else {
                        printf("    ");
                }

                /* Print the cell's score */
                if (optimal_path == 1) {
                        set_fmt(opt_path_fmt);
                }
                printf("%+*d", col_width, S->cells[col][row].score);
                if (optimal_path == 1) {
                        reset_fmt();
                }
        }
        printf("\n");
}

static void
print_table_row(score_table_t *S,
                walk_table_t *W,
                int row,
                char *s1,
                char *s2,
                int col_width,
                int unicode)
{
        print_directional_row(W, row, s1, s2, col_width, unicode);
        print_score_row(S, W, row, s1, s2, col_width, unicode);
}

static void
print_top_string(score_table_t *S, char *s1, int col_width)
{
        set_fmt(top_string_fmt);
        printf("*    %*s", col_width, "-");
        for (int i = 0; i < S->M - 1; i++) {
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
        return w + 1; /* add 1 to make room for a positive/negative sign */
}

void
print_table(score_table_t *S, walk_table_t *W, char *s1, char *s2, int unicode)
{
        int col_width = width_needed_to_print_integer(S->greatest_abs_val);

        /* Print the top string, i.e. the first input string. */
        print_top_string(S, s1, col_width);

        /* Print the rest of the rows, bordered on the left by the side
           string, i.e. the second input string. */
        for (int i = 0; i < S->N; i++) {
                print_table_row(S, W, i, s1, s2, col_width, unicode);
        }
}
