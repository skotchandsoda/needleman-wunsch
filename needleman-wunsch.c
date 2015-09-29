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
 * needleman-wunsch.c - Implementation of the Needleman-Wunsch sequence
 *                      alignment algorithm.
 *                      http://en.wikipedia.org/Needleman–Wunsch_algorithm
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
#include "read-sequences.h"
#include "table.h"

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
usage: needleman-wunsch [-c][-h][-l][-q][-s][-t][-u]\n\
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

void
print_aligned_strings_and_counts(char *X,
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
                if (no_print_strings == 0) {
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

        if (no_print_strings == 0) {
                printf("\n");

                for (int i = n; i > -1; i--) {
                        print_aligned_string_char(Y, X, i);
                }
                printf("\n");
        }

        // Print match/mismatch/gap counts if lflag was set
        if (print_counts == 1) {
                printf("%d match%s, %d mismatch%s, %d indel%s\n",
                       match_count, (match_count == 1 ? "" : "es"),
                       mismatch_count, (mismatch_count == 1 ? "" : "es"),
                       gap_count, (gap_count == 1 ? "" : "s"));
        }

        printf("\n");
}

void
inc_solution_count(computation_t *C)
{
        if (C->num_threads > 1) {
                pthread_rwlock_wrlock(&C->solution_count_rwlock);
        }

        C->solution_count = C->solution_count + 1;

        if (C->num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }
}

unsigned int
get_solution_count(computation_t *C)
{
        if (C->num_threads > 1) {
                pthread_rwlock_rdlock(&C->solution_count_rwlock);
        }

        unsigned int count = C->solution_count;

        if (C->num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }

        return count;
}

/* Do the walk iteratively because we'll overrun the stack on a
   sufficiently large input.  Yes, it is ugly, but it is necessary if we
   want to handle arbitrarily large inputs. */
void
walk_table(computation_t *C, char *X, char *Y, int start_i, int start_j)
{
        /* We move through the table starting at the bottom-right corner */
        table_t *T = C->scores_table;
        int i = T->M - 1;  // position (x direction)
        int j = T->N - 1;  // position (y direction)
        int n = 0;         // character count

        int rightmost_col = i;
        int bottommost_row = j;

        debug("Starting scores table walk.");
        while (!(i == rightmost_col &&
                 j == bottommost_row &&
                 1 == T->cells[i][j].up_done &&
                 1 == T->cells[i][j].diag_done &&
                 1 == T->cells[i][j].left_done)) {

                /* We've visited the cell, so mark it as part of the
                   optimal path */
                if (tflag == 1) {
                        T->cells[i][j].in_optimal_path = 1;
                }

                /* Special Case: We've reached the top-left corner of
                   the table.  Print the current solution. */
                if (i == 0 && j == 0) {
                        if (qflag != 1 || lflag == 1) {
                                print_aligned_strings_and_counts(X, Y, n-1,
                                                                 qflag, lflag);
                        }
                        inc_solution_count(C);
                }

                /* Base Case: Done in current cell.  Return to source
                 * cell. */
                if (T->cells[i][j].up_done &&
                    T->cells[i][j].diag_done &&
                    T->cells[i][j].left_done) {
                        /* Mark all possible paths as "not done" for
                           future visits */
                        T->cells[i][j].up_done = (T->cells[i][j].up ? 0 : 1);
                        T->cells[i][j].diag_done = (T->cells[i][j].diag ? 0 : 1);
                        T->cells[i][j].left_done = (T->cells[i][j].left ? 0 : 1);

                        /* Change i and j so we are "back in the source
                           cell."  Mark the source cell's relevant
                           direction "done" */
                        switch(T->cells[i][j].src_direction) {
                        case up:
                                j = j + 1;
                                T->cells[i][j].up_done = 1;
                                break;
                        case left:
                                i = i + 1;
                                T->cells[i][j].left_done = 1;
                                break;
                        case diag:
                                i = i + 1;
                                j = j + 1;
                                T->cells[i][j].diag_done = 1;
                                break;
                        default:
                                unreachable();
                        }
                        /* Decrement n so we can write in another
                           equivalent solution in a later pass */
                        n = n - 1;
                }
                /* Recursive Case: Not done in current cell.  Do stuff
                   in adjacent (up/diag/left) cells. */
                else {
                        if (1 == T->cells[i][j].diag &&
                            0 == T->cells[i][j].diag_done) {
                                X[n] = C->top_string[i-1];
                                Y[n] = C->side_string[j-1];
                                i = i-1;
                                j = j-1;
                                T->cells[i][j].src_direction = diag;
                        } else if (1 == T->cells[i][j].left &&
                                   0 == T->cells[i][j].left_done) {
                                X[n] = C->top_string[i-1];
                                Y[n] = '-';
                                i = i-1;
                                T->cells[i][j].src_direction = left;
                        } else if (1 == C->scores_table->cells[i][j].up &&
                                   0 == T->cells[i][j].up_done) {
                                X[n] = '-';
                                Y[n] = C->side_string[j-1];
                                j = j-1;
                                T->cells[i][j].src_direction = up;
                        }

                        n = n+1;
                }
        }

        debug("Finished scores table walk.");
}

struct walk_table_args {
        computation_t *C;
        char *X;
        char *Y;
        int start_i;
        int start_j;
};

void
mark_optimal_path_in_table(computation_t *C)
{
        int max_aligned_strlen;
        char *X;
        char *Y;

        /* Allocate buffers for printing the optimally aligned strings.  In the
           worst case they will need to be M+N characters long. */
        max_aligned_strlen = C->scores_table->M + C->scores_table->N;
        X = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        check(NULL != X, "malloc failed");
        Y = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        check(NULL != Y, "malloc failed");

        debug("Allocated temporary solution printing strings X and Y.");

        /* We walk through the table starting at the bottom-right-hand corner */
        int i = C->scores_table->M - 1;  // position (x direction)
        int j = C->scores_table->N - 1;  // position (y direction)

        /* Walk the table starting at the bottom-right corner, marking cells in
           the optimal path and counting the total possible optimal solutions
           (alignments) */
        walk_table(C, X, Y, i, j);

        /* Clean up buffers */
        free(X);
        free(Y);
}

void
inc_branch_count(computation_t *C)
{
        if (C->num_threads > 1) {
                pthread_rwlock_wrlock(&C->scores_table->branch_count_rwlock);
        }

        C->scores_table->branch_count = C->scores_table->branch_count + 1;

        if (C->num_threads > 1) {
                pthread_rwlock_unlock(&C->scores_table->branch_count_rwlock);
        }
}

static int
max3(int a, int b, int c)
{
        int m = a;
        (m < b) && (m = b);
        (m < c) && (m = c);
        return m;
}

void
process_cell(computation_t *C, int col, int row)
{
        /* Cell we want to compute the score for */
        cell_t *target_cell = &C->scores_table->cells[col][row];

        /* Cells we'll use to compute target_cell's score */
        cell_t *up_cell   = &C->scores_table->cells[col][row-1];
        cell_t *diag_cell = &C->scores_table->cells[col-1][row-1];
        cell_t *left_cell = &C->scores_table->cells[col-1][row];

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

        /*
         * BEGIN CRITICAL SECTIONS
         */

        if (C->num_threads > 1) {
                /* Wait for signal that left_cell is processed, then
                   lock the left cell's score mutex. */
                pthread_mutex_lock(&left_cell->score_mutex);
                while (left_cell->processed != 1) {
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

        /* The current cell's score is the max of the three candidate scores */
        target_cell->score = max3(up_score, left_score, diag_score);

        if (C->num_threads > 1) {
                /* We've finalized the cell's score, so we mark it as
                 * processed for any waiting thread. */
                target_cell->processed = 1;

                /* Signal the waiting thread and release the score mutex */
                pthread_cond_signal(&target_cell->processed_cv);
                pthread_mutex_unlock(&target_cell->score_mutex);
        }

        /*
         * END CRITICAL SECTIONS
         */

        /* Mark the optimal paths.  Provided that a path's score is
           equal to the target cell's score, i.e. the maximum of the
           three candidate scores, it is an optimal path. */
        if (target_cell->score == diag_score) {
                target_cell->diag = 1;
                target_cell->diag_done = 0;
        } else {
                target_cell->diag_done = 1;
        }
        if (target_cell->score == up_score) {
                target_cell->up = 1;
                target_cell->up_done = 0;
        } else {
                target_cell->up_done = 1;
        }
        if (target_cell->score == left_score) {
                target_cell->left = 1;
                target_cell->left_done = 0;
        } else {
                target_cell->left_done = 1;
        }

        /* If we can branch here, i.e. multiple paths have the same
           scores, note it. */
        if (target_cell->diag + target_cell->up + target_cell->left > 1) {
                inc_branch_count(C);
        }
}

void
process_column(computation_t *C, int col)
{
        table_t *T = C->scores_table;

        /* Compute the score for each cell in the column */
        for (int row = 1; row < C->scores_table->N; row++) {
                // Compute the cell's score
                process_cell(C, col, row);

                // If we're printing the table and the absolute value of
                // the current cell's score is greater than the one
                // marked in the table, update the largest value
                int current_abs_score = abs(T->cells[col][row].score);
                if (tflag == 1 && current_abs_score > T->greatest_abs_val) {
                        T->greatest_abs_val = current_abs_score;
                }
        }
}

void *
process_column_set(void *args)
{
        /* Unpack arguments.  We pass them in a struct because pthreads
           only lets us pass a block of memory as argument to the
           initial function */
        struct process_col_set_args *A = (struct process_col_set_args *)args;
        int current_col = A->start_col;
        computation_t *C = A->C;

        /* Process all columns in the thread's column set */
        while (current_col < C->scores_table->M) {
                process_column(C, current_col);
                current_col = current_col + C->num_threads;
        }

        return NULL; /* FIXME: Return some value indicating success? */
}

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
        for (int i = 0; i < C->num_threads; i++) {
                /* Initialize thread-local arguments for processing a
                   column set */
                args[i].start_col = i + 1;
                args[i].C = C;

                /* Spawn the thread */
                int res = pthread_create(&(C->worker_threads[i]), NULL,
                                         process_column_set,
                                         &args[i]);
                check(0 == res, "pthread_create failed");
        }

        /* Join the worker threads */
        int res;
        for (int i = 0; i < C->num_threads; i++) {
                res = pthread_join(C->worker_threads[i], NULL);
                check(0 == res, "pthread_join failed");
                debug("Joined thread %d", i+1);
        }
        debug("Joined %d worker thread%s", C->num_threads,
              (C->num_threads == 1 ? "" : "s"));
        free(C->worker_threads);
        debug("%u branches in table\n", C->scores_table->branch_count);
}



/* Print details about the algorithm's run, i.e. number of optimal
   alignments and maximum possible score */
void
print_summary(computation_t *C)
{
        unsigned int soln_count = get_solution_count(C);
        int max_col = C->scores_table->M - 1;
        int max_row = C->scores_table->N - 1;
        printf("%d optimal alignment%s\n",
               soln_count, (soln_count > 1 ? "s" : ""));
        printf("Optimal score is %-d\n",
               C->scores_table->cells[max_col][max_row].score);
}

void
needleman_wunsch(char *s1, char *s2, int m, int k, int d, int num_threads)
{
        /* Allocate and initialize computation */
        computation_t *C = alloc_computation();
        init_computation(C, s1, s2, m, k, d, num_threads);

        /* Fill out table, i.e. compute the optimal score */
        compute_table_scores(C);

        /* Walk the table.  Mark the optimal path if tflag is set, print
           the aligned strings if qflag is NOT set, and list counts for
           each alignment if lflag is set */
        if (qflag != 1 || lflag == 1 || sflag == 1 || tflag == 1) {
                mark_optimal_path_in_table(C);
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
                print_table(C->scores_table, C->top_string, C->side_string, uflag);
        }

        /* Clean up */
        free_computation(C);
}

int
main(int argc, char **argv)
{
        /* Strings to align */
        char *s1;
        char *s2;

        char *infile_path = NULL;
        FILE *in = NULL;

        /* Scoring values */
        int m, k, d;

        int num_threads = 1;

        // Parse option flags
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

        /* Set program name */
        set_prog_name(argv[0]);

        /* Make sure we have the right number of operands */
        if (optind + NUM_OPERANDS != argc) {
                log_err("expected %d operands; received %s %d", NUM_OPERANDS,
                        (argc - optind > NUM_OPERANDS ||
                         argc - optind == 0 ? "" : "only"),
                        argc - optind);
                usage();
        }

        /* If we got a filename, read the strings from that file.
           Otherwise, read the strings from stdin. */
        if (NULL == infile_path) {
                in = stdin;
        } else {
                in = fopen(infile_path, "r");
                check(NULL != in, "Failed to open %s", infile_path);
        }

        read_two_sequences_from_stream(&s1, &s2, in);

        /* Scoring values */
        m = atoi(argv[optind + 0]);
        k = atoi(argv[optind + 1]);
        d = atoi(argv[optind + 2]);

        /* Solve */
        needleman_wunsch(s1, s2, m, k, d, num_threads);

        /* Clean up */
        free(s1);
        free(s2);

        return 0;
}
