/*-
 * Copyright (c) 2015, Scott Cheloha.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following
 *      disclaimer in the documentation and/or other materials provided
 *      with the distribution.
 *
 *   3. Neither the name of the copyright holder nor the names of its
 *      contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* 
 * smith-waterman.c - an implementation of the Smith-Waterman algorithm
 *                    for locally optimal sequence alignment.
 *                    https://en.wikipedia.org/wiki/Smith–Waterman_algorithm
 */

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "computation.h"
#include "dbg.h"
#include "format.h"
#include "needleman-wunsch.h"
#include "print-table.h"
#include "read-sequences.h"
#include "score-table.h"
#include "walk-table.h"

#define GAP_CHAR '-'

#define NUM_OPERANDS 3

/* ANSI terminal output formatting flag (defined in format.h) */
extern int cflag;

/* Program name (defined in dbg.h) */
extern char *prog;

void
usage()
{
        fprintf(stderr, "\
usage: smith-waterman [-c][-h][-l][-q][-s][-t][-u]\n\
                      [-p num-threads] [-f sequence-file] m k d\n\
Align two sequences with the Needleman-Wunsch algorithm\n\
operands:\n\
   m   match bonus\n\
   k   mismatch penalty\n\
   d   indel (gap) penalty\n\
options:\n\
  -c   color the output with ANSI escape sequences\n\
  -f sequence-file\n\
       read the input strings from 'sequence-file' instead of standard input\n\
  -h   print this usage message\n\
  -l   list match, mismatch, and indel counts for each alignment pair\n\
  -p num-threads\n\
       parallelize the computation with 'num-threads' threads (must be >1)\n\
  -q   be quiet and don't print the aligned strings\n\
  -s   summarize the algorithm's run\n\
  -t   print the scores table; only useful for shorter input strings\n\
  -u   use unicode arrows when printing the scores table\n");
        exit(1);
}

/*
 * print_aligned_string_char()
 *
 *   Print the character s1[n] formatted according to its relationship
 *   to s2[n].  Depending on whether the two match, mismatch, or are gap
 *   characters, set the output formatting accordingly.
 *
 *   See format.h for definitions of the formats referenced in the
 *   function.
 */
void
print_aligned_string_char(char *s1, char *s2, int n)
{
        /* Format the output character as defined in format.h */
        if (s1[n] == s2[n]) {
                set_fmt(match_char_fmt);
        } else if (s1[n] == GAP_CHAR || s2[n] == GAP_CHAR) {
                set_fmt(gap_char_fmt);
        } else if (s1[n] != s2[n]) {
                set_fmt(mismatch_char_fmt);
        } else {
                unreachable();
        }

        /* Print the character */
        printf("%c", s1[n]);

        reset_fmt();
}

/*
 * print_aligned_strings_and_counts()
 *
 *   Print the aligned sequences X and Y unless no_print_strings is 1.
 *
 *   X - Aligned form of the top string to print
 *
 *   Y - Aligned form of the side string to print
 *
 *   n - Length in characters of X and Y
 *
 *   no_print_strings - If equal to 1, we don't print X and Y
 *
 *   print_counts - If equal to 1, we print match/mismatch/indel counts
 *                  for this pair of aligned sequences
 */
void
print_solution(char *X,
               char *Y,
               int n,
               int no_print_strings,
               int print_counts)
{
        int match_count = 0;
        int mismatch_count = 0;
        int gap_count = 0;

        /* Print the strings backwards */
        for (int i = n; i > -1; i--) {
                if (no_print_strings != 1) {
                        print_aligned_string_char(X, Y, i);
                }
                if (print_counts == 1) {
                        if (X[i] == Y[i]) {
                                match_count = match_count + 1;
                        } else if (X[i] == GAP_CHAR || Y[i] == GAP_CHAR) {
                                gap_count = gap_count + 1;
                        } else {
                                mismatch_count = mismatch_count + 1;
                        }
                }
        }

        if (0 == no_print_strings) {
                printf("\n");

                for (int i = n; i > -1; i--) {
                        print_aligned_string_char(Y, X, i);
                }
                printf("\n");
        }

        /* Print match/mismatch/gap counts if lflag was set */
        if (print_counts == 1) {
                printf("%d match%s, %d mismatch%s, %d indel%s\n",
                       match_count, (match_count == 1 ? "" : "es"),
                       mismatch_count, (mismatch_count == 1 ? "" : "es"),
                       gap_count, (gap_count == 1 ? "" : "s"));
        }

        printf("\n");
}

int
fill_rest_of_solution_buffers(computation_t *C,
                              char *X,
                              char *Y,
                              int i,
                              int j,
                              int n)
{
        /* While we aren't yet at the edge of the table, keep filling
         * the buffers with characters from the top/side strings */
        while (0 < i && 0 < j) {
                X[n] = C->top_string[i-1];
                Y[n] = C->side_string[j-1];
                n = n + 1;
                i = i - 1;
                j = j - 1;
        }

        /* If we're now at the top of the table, copy the rest of the
         * top string into buffer X and fill the remainder of Y with
         * spaces. */
        while (0 < i) {
                X[n] = C->top_string[i-1];
                Y[n] = ' ';
                n = n + 1;
                i = i - 1;
        }

        /* If we're now at the side of the table, copy the rest of the
         * side string into buffer Y and fill the remainder of X with
         * spaces. */
        while (0 < j) {
                X[n] = ' ';
                Y[n] = C->side_string[j-1];
                n = n + 1;
                j = j - 1;
        }

        return n;
}

/*
 * construct_alignments_from_cell()
 *
 *   Starting at cell (start_i, start_j), iterate through the given
 *   computation's walk_table and reconstruct all optimal alignments of
 *   the input strings.  The cell (start_i, start_j) forms the
 *   bottom-righthand boundary of the subtable this call will construct
 *   solutions for.
 *
 *   C - computation instance to reconstruct alignments for
 *
 *   X - buffer to store the aligned top string in
 *
 *   Y - buffer to store the aligned side string in
 *
 *   start_col - column of the cell to begin iterating from, i.e. the
 *               right boundary column of the subtable we're
 *               constructing solutions for
 *
 *   start_row - row of the cell to begin iterating from, i.e. the lower
 *               boundary row of the subtable we're constructing
 *               solutions for
 *
 *   start_char_count - starting offset in the alignment string buffers
 *                      (X & Y)
 */
void
construct_alignments_for_subtable(computation_t *C,
                                  char *X,
                                  char *Y,
                                  int start_col,
                                  int start_row,
                                  int start_char_count)
{
        /* We move through the walk table starting at the bottom-right
         * corner as defined by start_i (the righthand limit for this
         * table walk) and start_j (the lower limit for this table
         * walk). */
        walk_table_t *W = C->walk_table;
        int i = start_col;  /* position (x direction) */
        int j = start_row;  /* position (y direction) */
        int n = start_char_count;  /* character count */

        debug("Starting alignment construction from (%d,%d)...", i, j);

        /* We do the walk iteratively because we'll overrun the stack on
         * a sufficiently large input.  Yes, it is ugly, but it is
         * necessary if we want to handle arbitrarily large inputs. */
        while (0 <= i && 0 <= j && start_col >= i && start_row >= j &&
               !(i == start_col &&
                 j == start_row &&
                 1 == W->cells[i][j].up_done &&
                 1 == W->cells[i][j].diag_done &&
                 1 == W->cells[i][j].left_done)) {

                /* We've visited the cell, so mark it as part of the
                 * optimal path */
                if (tflag == 1) {
                        W->cells[i][j].in_optimal_path = 1;
                }

                /*
                 * Special Case: We can go no further from the current
                 *               cell, i.e.  no forward paths are
                 *               possible, so we print the current
                 *               solution (the aligned strings X & Y) to
                 *               the standard output.
                 */
                if (W->cells[i][j].up == 0 &&
                    W->cells[i][j].diag == 0 &&
                    W->cells[i][j].left == 0) {
                        /* If we haven't exhausted the top/side strings,
                         * fill out the rest of X and Y. */
                        int new_n = fill_rest_of_solution_buffers(C, X, Y, i, j, n);

                        /* Then, print the solution and increment the
                         * solution total. */
                        if (qflag != 1 || lflag == 1) {
                                print_solution(X, Y, new_n-1, qflag, lflag);
                        }

                        inc_solution_count(C);
                }

                /*
                 * Base Case: All cells adjacent (up/diag/left) to the
                 *            current cell have been marked "done," so
                 *            we (1) print the current solution, and
                 *            then (2) return to the cell we were last
                 *            in via the 'src_direction' indicator.
                 */
                if (W->cells[i][j].up_done &&
                    W->cells[i][j].diag_done &&
                    W->cells[i][j].left_done) {
                        /* Mark all possible paths as "not done" for
                           future visits. */
                        W->cells[i][j].up_done   = (W->cells[i][j].up ? 0 : 1);
                        W->cells[i][j].diag_done = (W->cells[i][j].diag ? 0 : 1);
                        W->cells[i][j].left_done = (W->cells[i][j].left ? 0 : 1);

                        /* Change i and j so we are "back in the source
                           cell."  Mark the source cell's relevant
                           direction "done." */
                        switch(W->cells[i][j].src_direction) {
                        case up:
                                j = j + 1;
                                W->cells[i][j].up_done = 1;
                                break;
                        case left:
                                i = i + 1;
                                W->cells[i][j].left_done = 1;
                                break;
                        case diag:
                                i = i + 1;
                                j = j + 1;
                                W->cells[i][j].diag_done = 1;
                                break;
                        default:
                                unreachable();
                        }

                        /* Decrement n so we can write in another
                         * equivalent solution in a later pass */
                        n = n - 1;
                }

                /*
                 * Recursive Case: Not done in current cell.  Copy
                 *                 characters from top/side strings into
                 *                 X and Y as needed and then iterate
                 *                 into an adjacent (up/diag/left) cell
                 *                 if we haven't yet marked the cell
                 *                 "done."
                 */
                else {
                        if (1 == W->cells[i][j].diag &&
                            0 == W->cells[i][j].diag_done) {
                                X[n] = C->top_string[i-1];
                                Y[n] = C->side_string[j-1];
                                i = i - 1;
                                j = j - 1;
                                W->cells[i][j].src_direction = diag;
                        } else if (1 == W->cells[i][j].left &&
                                   0 == W->cells[i][j].left_done) {
                                X[n] = C->top_string[i-1];
                                Y[n] = '-';
                                i = i - 1;
                                W->cells[i][j].src_direction = left;
                        } else if (1 == W->cells[i][j].up &&
                                   0 == W->cells[i][j].up_done) {
                                X[n] = '-';
                                Y[n] = C->side_string[j-1];
                                j = j - 1;
                                W->cells[i][j].src_direction = up;
                        }

                        n = n + 1;
                }
        }

        debug("Finished alignment construction from (%d,%d).", i, j);
}

struct row_col_node {
        int col;
        int row;
        struct row_col_node *next;
};

struct row_col_node *
alloc_row_col_node()
{
        struct row_col_node *n =
                (struct row_col_node *)malloc(sizeof(struct row_col_node));
        check(NULL != n, "malloc failed");
        n->next = NULL;

        return n;
}

struct row_col_node_list {
        int len;
        struct row_col_node *head;
};

struct row_col_node_list *
alloc_row_col_node_list()
{
        struct row_col_node_list *L =
                (struct row_col_node_list *)malloc(sizeof(struct row_col_node_list));
        check(NULL != L, "malloc failed");
        L->len = 0;
        L->head = NULL;

        return L;
}

void
free_row_col_node_list(struct row_col_node_list *L)
{
        struct row_col_node *curr = L->head;
        struct row_col_node *prev = curr;
        while (NULL != curr) {
                prev = curr;
                curr = curr->next;
                free(prev);
        }

        free(L);
}

struct row_col_node_list *
get_list_of_starting_cells(computation_t *C)
{
        /* Allocate list */
        struct row_col_node_list *L = alloc_row_col_node_list();
        struct row_col_node **ptr_to_end_of_L = &L->head;

        /* Populate list */
        struct row_col_node *node;
        score_table_t *S = C->score_table;
        int max_score = S->greatest_abs_val;
        debug("Looking for start-cells for optimal local alignments.  " \
              "Max score is %d.", max_score);
        for (int i = 1; i < S->M; i++) {
                for (int j = 1; j < S->N; j++) {
                        if (max_score == S->cells[i][j].score) {
                                debug("Start cell for local alignment " \
                                      "@ (%d,%d)", i, j);
                                node = alloc_row_col_node();
                                node->col = i;
                                node->row = j;
                                *ptr_to_end_of_L = node;
                                ptr_to_end_of_L = &(*ptr_to_end_of_L)->next;
                                L->len = L->len + 1;
                        }
                }
        }
        debug("Found %d eligible cells.", L->len);

        return L;
}

/*
 * construct_alignments()
 *
 *   Construct all optimal alignments for the walk table of the given
 *   computation instance.  It the '-q' flag is not set, all optimal
 *   alignments will be printed to the standard output.
 *
 *   NOTE: This routine could be modified to run
 *         construct_alignments_for_subtable() in parallel, which would
 *         significantly improve performance for large/dissimilar input
 *         strings (i.e. strings that produce computations with hundreds
 *         of thousands of branches in their reference walk table).
 *
 *   C - computation instance to construct optimal alignments for
 */
void
construct_alignments(computation_t *C)
{
        struct row_col_node_list *L = get_list_of_starting_cells(C);

        int max_aligned_strlen;
        char *X;
        char *Y;

        /* Allocate buffers for printing the optimally aligned strings.
         * In the worst case they will need to be M+N characters
         * long. */
        max_aligned_strlen = C->score_table->M + C->score_table->N;

        debug("Allocating temporary solution printing strings X and Y.");
        X = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        check(NULL != X, "malloc failed");
        Y = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        check(NULL != Y, "malloc failed");

        /* For each node in the list, we walk through the table starting
         * at the cell indicated by the node */
        struct row_col_node *node = L->head;
        while (NULL != node) {
                int start_col = node->col;
                int start_row = node->row;
                int n = 0;

                int i = C->score_table->M - 1;
                int j = C->score_table->N - 1;

                /* Adjust horizontally on walk table */
                while (i - start_col > j - start_row) {
                        X[n] = C->top_string[i-1];
                        Y[n] = ' ';
                        i = i - 1;
                        n = n + 1;
                }

                /* Adjust vertically on walk table */
                while (j - start_row > i - start_col) {
                        X[n] = ' ';
                        Y[n] = C->side_string[j-1];
                        j = j - 1;
                        n = n + 1;
                }

                /* Move diagonally on walk table until at starting col/row */
                while (j != start_row || i != start_col) {
                        X[n] = C->top_string[i-1];
                        Y[n] = C->side_string[j-1];
                        i = i - 1;
                        j = j - 1;
                        n = n + 1;
                }

                /* Walk the cell at (start_col, start_row), marking
                 * cells in the optimal path and counting the total
                 * possible optimal solutions (alignments) */
                construct_alignments_for_subtable(C, X, Y, start_col, start_row, n);

                node = node->next;
        }

        /* Clean up solution storage buffers */
        free(X);
        free(Y);

        free_row_col_node_list(L);
}

/*
 * max3()
 *
 *   Return the maximum of the set {a, b, c}.
 */
static int
max3(int a, int b, int c)
{
        int m = a;
        if (m < b)
                m = b;
        if (m < c)
                m = c;
        return m;
}

/*
 * max3_or_zero()
 *
 *   Return the maximum of the set {a, b, c, 0}.
 */
static int
max3_or_zero(int a, int b, int c)
{
        int m = a;
        if (m < b)
                m = b;
        if (m < c)
                m = c;
        if (m < 0)
                m = 0;
        return m;
}

/*
 * score_cell()
 *
 *   Write the alignment score to the score_table cell at (col,row).
 *
 *   C - pointer to computation_t instance containing the target
 *       score_table
 *
 *   col - column of the target cell in the score table
 *
 *   row - row of the target cell in the score table
 */
void
score_cell(computation_t *C, int col, int row)
{
        /* Cell we want to compute the score for */
        score_table_cell_t *target_cell = &C->score_table->cells[col][row];

        /* Cells we'll use to compute target_cell's score */
        score_table_cell_t *up_cell   = &C->score_table->cells[col][row-1];
        score_table_cell_t *diag_cell = &C->score_table->cells[col-1][row-1];
        score_table_cell_t *left_cell = &C->score_table->cells[col-1][row];

        /* Candidate scores */
        int up_score = up_cell->score - C->indel_penalty;
        int diag_score = 0;
        if (C->top_string[col-1] == C->side_string[row-1]) {
                diag_score = diag_cell->score + C->match_score;
                target_cell->match = 1;
        } else {
                diag_score = diag_cell->score - C->mismatch_penalty;
                target_cell->match = 0;
        }

        /***************************
         * BEGIN CRITICAL SECTIONS *
         ***************************/

        if (C->num_threads > 1) {
                /* Wait for signal that left_cell is processed, then
                   lock the left cell's score mutex. */
                pthread_mutex_lock(&left_cell->score_mutex);
                while (1 != left_cell->processed) {
                        pthread_cond_wait(&left_cell->processed_cv,
                                          &left_cell->score_mutex);
                }
        }

        int left_score = left_cell->score - C->indel_penalty;

        if (C->num_threads > 1) {
                /* We're done with the left score, so free the mutex */
                pthread_mutex_unlock(&left_cell->score_mutex);

                /* Lock current cell's score mutex and process the cell */
                pthread_mutex_lock(&target_cell->score_mutex);
        }

        /* For Needleman-Wunsch and Overlap, the current cell's score is
         * the max of the three candidate scores.  For Smith-Waterman,
         * the current cell's score is the max of the aforementioned
         * cell-scores and zero. */
        if (C->algorithm == SW) {
                target_cell->score = max3_or_zero(up_score, left_score, diag_score);
        } else {
                target_cell->score = max3(up_score, left_score, diag_score);
        }

        if (C->num_threads > 1) {
                /* We've finalized the cell's score, so we mark it as
                 * processed for any waiting thread. */
                target_cell->processed = 1;

                /* Signal the waiting thread and release the score mutex */
                pthread_cond_signal(&target_cell->processed_cv);
                pthread_mutex_unlock(&target_cell->score_mutex);
        }

        /*************************
         * END CRITICAL SECTIONS *
         *************************/

        /* Mark the optimal paths in the walk table.  Provided that a
         * path's score is equal to the target cell's score, i.e. the
         * maximum of the three candidate scores, it is an optimal
         * path. */
        walk_table_cell_t *target_walk_cell = &C->walk_table->cells[col][row];

        /* If we're using Smith-Waterman and the current cell's score is zero,
         * this is the end of the "path." */
        if (C->algorithm == SW && target_cell->score == 0) {
                target_walk_cell->diag_done = 1;
                target_walk_cell->up_done = 1;
                target_walk_cell->left_done = 1;
        }
        /* Otherwise, mark the path accordingly. */
        else {
                if (target_cell->score == diag_score) {
                        target_walk_cell->diag = 1;
                        target_walk_cell->diag_done = 0;
                } else {
                        target_walk_cell->diag_done = 1;
                }
                if (target_cell->score == up_score) {
                        target_walk_cell->up = 1;
                        target_walk_cell->up_done = 0;
                } else {
                        target_walk_cell->up_done = 1;
                }
                if (target_cell->score == left_score) {
                        target_walk_cell->left = 1;
                        target_walk_cell->left_done = 0;
                } else {
                        target_walk_cell->left_done = 1;
                }
        }
        /* If we can branch here, i.e. multiple paths have the same
           score, note it by incrementing the branch count. */
        if (1 < (target_walk_cell->diag +
                 target_walk_cell->up +
                 target_walk_cell->left)) {
                inc_branch_count(C->walk_table, C->num_threads);
        }
}

/*
 * score_cell_column()
 *
 *   Write alignment scores to a column of cells in a computation's
 *   score table.
 *
 *     C - pointer to the computation instance containing the target
 *         score table
 *
 *   col - index of the column of cells to score
 */
void
score_cell_column(computation_t *C, int col)
{
        score_table_t *S = C->score_table;

        /* Compute the score for each cell in the column */
        for (int row = 1; row < S->N; row++) {
                /* Compute the cell's score */
                score_cell(C, col, row);

                /*
                 * If we're printing the table and the absolute value of
                 * the current cell's score is greater than the one
                 * marked in the table, update the largest value.
                 */
                int current_abs_score = abs(S->cells[col][row].score);
                if (tflag == 1 && current_abs_score > S->greatest_abs_val) {
                        S->greatest_abs_val = current_abs_score;
                }
        }
}

/*
 * score_cell_column_set()
 *
 *   Write alignment scores to a set of cell columns in a computation's
 *   score table.  Given a starting column x, the current thread will score
 *   columns x + i*num_threads for i=0 until x + i*num_threads exceeds the
 *   total number of columns in the table.
 *
 *   args - pointer to a struct process_col_set_args, which contains a
 *          pointer to the target computation instance and a column
 *          index for the thread to start with
 */
void *
score_cell_column_set(void *args)
{
        /* Unpack arguments.  We pass them in a struct because pthreads
           only lets us pass a block of memory as argument to the
           initial function */
        struct process_col_set_args *A = (struct process_col_set_args *)args;
        int current_col = A->start_col;
        computation_t *C = A->C;

        /* Process all columns in the thread's column set */
        while (current_col < C->score_table->M) {
                score_cell_column(C, current_col);
                current_col = current_col + C->num_threads;
        }

        return NULL; /* FIXME: Return some value indicating success? */
}

/*
 * compute_table_scores()
 *
 *   Score each cell in a computation instance's score table.
 *
 *   C - target computation instance
 */
void
compute_table_scores(computation_t *C)
{
        /* Allocate storage for thread ids and arguments to process_col_set */
        C->worker_threads = (pthread_t *)malloc(C->num_threads * sizeof(pthread_t));
        check(NULL != C->worker_threads, "malloc failed");
        struct process_col_set_args *args;
        args = (struct process_col_set_args *)malloc(C->num_threads *
                                                     sizeof(struct process_col_set_args));
        check(NULL != args, "malloc failed");

        /* Spawn worker threads to process sets of columns */
        debug("Spawning %d worker thread%s for scores table computation",
              C->num_threads, (C->num_threads == 1 ? "" : "s"));
        for (unsigned int i = 0; i < C->num_threads; i++) {
                /* Initialize thread-local arguments for processing a
                 * set of cell-columns */
                args[i].start_col = i + 1;
                args[i].C = C;

                /* Spawn the thread */
                int res = pthread_create(&C->worker_threads[i],
                                         NULL,
                                         score_cell_column_set,
                                         &args[i]);
                check(0 == res, "pthread_create failed");
        }

        /* Join the worker threads */
        int res;
        unsigned int join_count = 0;
        for (unsigned int i = 0; i < C->num_threads; i++) {
                res = pthread_join(C->worker_threads[i], NULL);
                check(0 == res, "pthread_join failed");
                join_count = join_count + 1;
                debug("Joined thread %d", i+1);
        }
        check(join_count == C->num_threads, "this should never happen");
        debug("Joined %d worker thread%s", C->num_threads,
              (C->num_threads == 1 ? "" : "s"));
        free(C->worker_threads);
        debug("%u branches in walk table\n",
              get_branch_count(C->walk_table, C->num_threads));
}

/*
 * smith_waterman()
 *
 *   Execute the Smith-Waterman locally-optimal sequence alignment
 *   algorithm for the given inputs.
 *
 *   s1 - Top string, i.e. the first input sequence
 *
 *   s2 - Side string, i.e. the second input sequence
 *
 *    m - Match bonus, i.e. the amount added to the diagonal cell's
 *        score when it is optimal for the two characters in a cell to
 *        be the same (which, depending on the value of m, is probably
 *        true).
 *
 *    k - Mismatch penalty, i.e. the amount subtracted from the diagonal
 *        cell's score when it is optimal for the two characters in a
 *        cell to not match.
 *
 *    d - Indel penalty, i.e. the amount subtracted from the upper or
 *        left cell's score when it is optimal to skip the character in
 *        the opposite sequence.
 *
 *    num_threads - the number of threads to execute in parallel when
 *                  scoring the computation's score table.
 */
void
smith_waterman(char *s1, char *s2, int m, int k, int d, int num_threads)
{
        /* Allocate and initialize computation */
        computation_t *C = alloc_computation();
        init_computation(C, SW, s1, s2, m, k, d, num_threads);

        /* Fill out table, i.e. compute the optimal score */
        compute_table_scores(C);

        /* Walk the table.  Mark the optimal path if tflag is set, print
         * the aligned strings if qflag is NOT set, and list counts for
         * each alignment if lflag is set */
        if (qflag != 1 || lflag == 1 || sflag == 1 || tflag == 1) {
                construct_alignments(C);
        }

        /* Print summary if sflag is set */
        if (sflag == 1) {
                print_summary(C);
        }

        /* Print table if tflag is set */
        if (tflag == 1) {
                /* Print an extra newline to separate the output
                 * sections */
                if (qflag != 1 || sflag == 1 || lflag == 1) {
                        printf("\n");
                }
                print_table(C->score_table, C->walk_table,
                            C->top_string, C->side_string, uflag);
        }

        /* Clean up */
        free_computation(C);
}

/*
 * main()
 *
 *   Parse option arguments, read in the two input strings (s1 and s2),
 *   and execute the Needleman-Wunsch algorithm for the aforementioned
 *   input strings and the operands m, k, and d.
 */
int main(int argc, char **argv)
{
        /* Strings to align */
        char *s1;
        char *s2;

        char *infile_path = NULL;
        FILE *in = NULL;

        /* Scoring values */
        int m, k, d;

        int num_threads = 1;

        /* Set program name */
        set_prog_name(argv[0]);

        /* Clear errno */
        errno = 0;

        /* Parse option flags */
        extern char *optarg;
        extern int optind;
        int c;

        while ((c = getopt(argc, argv, "cf:hlp:qstu")) != -1) {
                switch (c) {
                case 'c':
                        cflag = 1;
                        break;
                case 'f':
                        infile_path = optarg;
                        break;
                case 'h':
                        usage();
                        break;
                case 'l':
                        lflag = 1;
                        break;
                case 'p':
                        num_threads = atoi(optarg);
                        check(num_threads > 1,
                              "num-threads == %d; num-threads "         \
                              "must be greater than 1", num_threads);
                        break;
                case 'q':
                        qflag = 1;
                        break;
                case 's':
                        sflag = 1;
                        break;
                case 't':
                        tflag = 1;
                        break;
                case 'u':
                        uflag = 1;
                        break;
                case '?':
                default:
                        usage();
                        break;
                }
        }

        /* Make sure we have the right number of operands */
        if (optind + NUM_OPERANDS != argc) {
                log_err("expected %d operands but received%s %d", NUM_OPERANDS,
                        (argc - optind > NUM_OPERANDS ||
                         argc - optind == 0 ? "" : " only"),
                        argc - optind);
                usage();
        }

        /* If we got a filename, read the strings from that file.
           Otherwise, read the strings from stdin. */
        if (NULL == infile_path) {
                in = stdin;
        } else {
                in = fopen(infile_path, "r");
                check(NULL != in, "failed to open %s", infile_path);
        }

        read_two_sequences_from_stream(&s1, &s2, in);

        /* Set scoring values to operands give on command-line */
        m = atoi(argv[optind + 0]);
        k = atoi(argv[optind + 1]);
        d = atoi(argv[optind + 2]);

        /* Solve the alignment */
        smith_waterman(s1, s2, m, k, d, num_threads);

        /* Clean up */
        free(s1);
        free(s2);

        return 0;
}
