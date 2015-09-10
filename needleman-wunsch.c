// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch - align two strings with the Needleman-Wunsch algorithm
//                    (see http://en.wikipedia.org/Needleman–Wunsch_algorithm)

#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "format.h"
#include "needleman-wunsch.h"
#include "table.h"

#define GAP_CHAR '-'

#define INPUT_STRING_BASE_SIZE 4096

/* Global flags defined in format.h */
extern int cflag;

void
usage()
{
        fprintf(stderr, "\
usage: needleman-wunsch [-c][-h][-l]\n\
           [-p num_threads][-q][-s][-t][-u] s1 s2 m k g\n\
Align two strings with the Needleman-Wunsch algorithm\n\
operands:\n\
  s1   first string to align (the top string)\n\
  s2   second string to align (the side string)\n\
   m   match bonus\n\
   k   mismatch penalty\n\
   g   gap penalty\n\
options:\n\
  -c   color the output with ANSI escape sequences\n\
  -h   print this usage message\n\
  -l   list match, mismatch, and gap counts after each alignment pair\n\
  -p num_threads\n\
       parallelize the computation with 'num_threads' threads (must be >1)\n\
  -q   be quiet and don't print the alignmened strings (cancels the '-l' flag)\n\
  -s   summarize the algorithm's run\n\
  -t   print the scores table; probably only useful for shorter input strings\n\
  -u   use unicode arrows when printing the scores table\n"
);
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
                fprintf(stderr, "the impossible has happened; giving up\n");
                exit(1);
        }

        /* Print the character */
        printf("%c", s1[n]);

        reset_fmt();
}

void
print_aligned_strings_and_counts(char *X, char *Y, int n, int print_counts)
{
        int match_count = 0;
        int mismatch_count = 0;
        int gap_count = 0;

        /* Print the strings backwards */
        for (int i = n; i > -1; i--) {
                print_aligned_string_char(X, Y, i);
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
        printf("\n");

        for (int i = n; i > -1; i--) {
                print_aligned_string_char(Y, X, i);
        }
        printf("\n");

        // Print match/mismatch/gap counts if iflag was set
        if (print_counts == 1) {
                printf("%d match%s, %d mismatch%s, %d gap%s\n",
                       match_count, (match_count == 1 ? "" : "es"),
                       mismatch_count, (mismatch_count == 1 ? "" : "es"),
                       gap_count, (gap_count == 1 ? "" : "s"));
        }
}

void
inc_solution_count(computation_t *C)
{
        if (num_threads > 1) {
                pthread_rwlock_wrlock(&C->solution_count_rwlock);
        }

        C->solution_count = C->solution_count + 1;

        if (num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }
}

unsigned int
get_solution_count(computation_t *C)
{
        if (num_threads > 1) {
                pthread_rwlock_rdlock(&C->solution_count_rwlock);
        }

        unsigned int count = C->solution_count;

        if (num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }

        return count;
}

/* Do the walk iteratively because we'll overrun the stack with new frames
   on a sufficiently large input.  Yes, it is ugly, but it is necessary
   if we want to handle arbitrarily large inputs. */
void
walk_table(computation_t *C, char *X, char *Y)
{
        // We move through the table starting at the bottom-right corner
        table_t *T = C->scores_table;
        int i = T->M - 1;  // position (x direction)
        int j = T->N - 1;  // position (y direction)
        int n = 0;         // character count

        int rightmost_col = i;
        int bottommost_row = j;

        fprintf(stderr, "starting walk\n");
        while (!(i == rightmost_col &&
                 j == bottommost_row &&
                 1 == T->cells[i][j].up_done &&
                 1 == T->cells[i][j].diag_done &&
                 1 == T->cells[i][j].left_done)) {
//                fprintf(stderr, "@ (%d,%d)\n", i, j);

                // We've visited the cell, so mark it as part of the
                // optimal path
                if (tflag == 1) {
                        T->cells[i][j].in_optimal_path = 1;
                }

                // Special Case: We've reached the top-left corner of the table.
                //               Print the current solution.
                if (i == 0 && j == 0) {
                        if (qflag != 1) {
                                if (get_solution_count(C) > 0) {
                                        printf("\n");
                                }
                                print_aligned_strings_and_counts(X, Y, n-1, lflag);
                        }
                        inc_solution_count(C);
                }

                // Base Case: Done in current cell.  Return to source cell.
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
                                fprintf(stderr, "the impossible has happened; " \
                                        "giving up\n");
                                exit(1);
                        }
                        /* Decrement n so we can write in another
                           equivalent solution in a later pass */
                        n = n - 1;
                }
                // Recursive Case: Not done in current cell.  Do stuff in
                //                 adjacent (up/diag/left) cells.
                else {
                        if (1 == T->cells[i][j].diag &&
                            0 == T->cells[i][j].diag_done) {
                                if (qflag != 1) {
                                        X[n] = C->top_string[i-1];
                                        Y[n] = C->side_string[j-1];
                                }
                                i = i-1;
                                j = j-1;
                                T->cells[i][j].src_direction = diag;
                                n = n+1;
                        } else if (1 == T->cells[i][j].left &&
                                   0 == T->cells[i][j].left_done) {
                                if (qflag != 1) {
                                        X[n] = C->top_string[i-1];
                                        Y[n] = '-';
                                }
                                i = i-1;
                                T->cells[i][j].src_direction = left;
                                n = n+1;
                        } else if (1 == C->scores_table->cells[i][j].up &&
                                   0 == T->cells[i][j].up_done) {
                                if (qflag != 1) {
                                        X[n] = '-';
                                        Y[n] = C->side_string[j-1];
                                }
                                j = j-1;
                                T->cells[i][j].src_direction = up;
                                n = n+1;
                        }
                }
        }
        fprintf(stderr, "done walk\n");
}

void
walk_table_recursively(computation_t *C,
                       char *X,
                       char *Y,
                       int i,
                       int j,
                       int n)
{
        // Mark the cell so that (later) the table printing
        // routine can mark the cell as part of an optimal path
        if (tflag == 1) {
                C->scores_table->cells[i][j].in_optimal_path = 1;
        }

        // Base case: We've reached the top-left corner of the table, i.e.
        //            the current cell_t has no direction bits set.
        if ((C->scores_table->cells[i][j].diag ||
             C->scores_table->cells[i][j].left ||
             C->scores_table->cells[i][j].up) == 0) {
                if (qflag != 1) {
                        if (get_solution_count(C) > 0) {
                                printf("\n");
                        }
                        print_aligned_strings_and_counts(X, Y, n-1, lflag);
                }
                inc_solution_count(C);
        }
        // Recursive case: We're not at the top-left corner of the table.
        //                 Make a recursive call for each direction bit set.
        else {
                fprintf(stderr, "@ (%d,%d)\n", i, j);
                if (C->scores_table->cells[i][j].diag == 1) {
                        if (qflag != 1) {
                                X[n] = C->top_string[i-1];
                                Y[n] = C->side_string[j-1];
                        }
                        walk_table_recursively(C, X, Y, i-1, j-1, n+1);
                }
                if (C->scores_table->cells[i][j].left == 1) {
                        if (qflag != 1) {
                                X[n] = C->top_string[i-1];
                                Y[n] = '-';
                        }
                        walk_table_recursively(C, X, Y, i-1, j, n+1);
                }
                if (C->scores_table->cells[i][j].up == 1) {
                        if (qflag != 1) {
                                X[n] = '-';
                                Y[n] = C->side_string[j-1];
                        }
                        walk_table_recursively(C, X, Y, i, j-1, n+1);
                }
        }
}


void
mark_optimal_path_in_table(computation_t *C)
{
        int max_aligned_strlen;
        char *X;
        char *Y;

        // If the qflag wasn't set, allocate buffers for printing the
        // optimally aligned strings.  In the worst case they will be
        // M+N characters long.
        if (qflag != 1) {
                max_aligned_strlen = C->scores_table->M + C->scores_table->N;
                X = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
                if (X == NULL) {
                        perror("malloc failed");
                        exit(1);
                }
                Y = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
                if (Y == NULL) {
                        perror("malloc failed");
                        exit(1);
                }
        }

        fprintf(stderr, "allocated temp printing strings X&Y\n");

        // We move through the table starting at the bottom-right corner
        int i = C->scores_table->M - 1;  // position (x direction)
        int j = C->scores_table->N - 1;  // position (y direction)
        int n = 0;                       // character count

        // Walk the table starting at the bottom-right corner recursively,
        // marking cells in the optimal path and counting the total possible
        // optimal solutions (alignments)
//        walk_table_recursively(C, X, Y, i, j, n);
        walk_table(C, X, Y);

        // Clean up buffers if we allocated for them
        if (qflag != 1) {
                free(X);
                free(Y);
        }
        // Print details about the algorithm's run, i.e. number of
        // optimal alignments and maximum possible score
        if (sflag == 1) {
                if (qflag != 1) {
                        printf("\n");
                }
                unsigned int soln_count = get_solution_count(C);
                printf("%d optimal alignment%s\n",
                       soln_count, (soln_count > 1 ? "s" : ""));
                printf("Optimal score is %-d\n",
                       C->scores_table->cells[i][j].score);
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
process_cell(table_t *T, int col, int row, char *s1, char *s2, int m, int k, int g)
{
        // Cell we want to compute the score for
        cell_t *target_cell = &T->cells[col][row];

        // Cells we'll use to compute target_cell's score
        cell_t *up_cell   = &T->cells[col][row-1];
        cell_t *diag_cell = &T->cells[col-1][row-1];
        cell_t *left_cell = &T->cells[col-1][row];

        // Candidate scores
        int up_score = up_cell->score - g;
        int diag_score = 0;
        if (s1[col-1] == s2[row-1]) {
                diag_score = diag_cell->score + m;
                target_cell->match = 1;
        } else {
                diag_score = diag_cell->score - k;
                target_cell->match = 0;
        }

        /*
         * BEGIN CRITICAL SECTIONS
         */

        if (num_threads > 1) {
                /* Lock current cell's score mutex and process the cell */
                pthread_mutex_lock(&target_cell->score_mutex);

                /* Wait for signal that left_cell is processed */
                pthread_mutex_lock(&left_cell->score_mutex);
                while (left_cell->processed == 0) {
                        pthread_cond_wait(&left_cell->processed_cv,
                                          &left_cell->score_mutex);
                }
        }

        int left_score = left_cell->score - g;

        if (num_threads > 1) {
                pthread_mutex_unlock(&left_cell->score_mutex);
        }

        // The current cell's score is the max of the three candidate scores
        target_cell->score = max3(up_score, left_score, diag_score);

        /* We've set the cell's score, so mark it as processed */
        target_cell->processed = 1;

        /* Signal the waiting thread and release the score mutex */
        if (num_threads > 1) {
                pthread_cond_signal(&target_cell->processed_cv);
                pthread_mutex_unlock(&target_cell->score_mutex);
        }

        /*
         * END CRITICAL SECTIONS
         */

        // Mark the optimal paths.  Provided that a path's
        // score is equal to the target cell's score, i.e. the maximum
        // of the three candidate scores, it is an optimal path.
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
}

void
process_column(table_t *T, int col, char *s1, char *s2, int m, int k, int g)
{
        fprintf(stderr, "processing col %d\n", col);
        // Compute the score for each cell in the column
        for (int row = 1; row < T->N; row++) {
                // Compute the cell's score
                process_cell(T, col, row, s1, s2, m, k, g);

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
process_column_set(void *start_col)
{
        int current_col = *((int *)start_col);
        while (current_col < comp->scores_table->M) {
                process_column(comp->scores_table, current_col,
                               comp->top_string, comp->side_string,
                               comp->match_score, comp->mismatch_penalty,
                               comp->gap_penalty);
                current_col = current_col + num_threads;
        }

        fprintf(stderr, "done processing col-set starting at col %d\n", *((int *)start_col));

        return NULL;
}

void
compute_table_scores()
{
        // If we're printing the table, initialize the largest value
        if (tflag == 1) {
                comp->scores_table->greatest_abs_val = 0;
        }

        // Allocate storage for thread ids and starting columns
        worker_threads = (worker_thread_t *)malloc(num_threads * sizeof(worker_thread_t));
        if (NULL == worker_threads) {
                fprintf(stderr, "malloc failed");
                exit(1);
        }

        // Spawn worker threads to process sets of columns
        for (int i = 0; i < num_threads; i++) {
                worker_threads[i].start_col = i + 1;
                int res = pthread_create(&worker_threads[i].thread_id, NULL,
                                         process_column_set,
                                         &worker_threads[i].start_col);
                if (0 != res) {
                        perror("pthread_create failed");
                        exit(1);
                }
        }

        fprintf(stderr, "spawned threads\n");

        // Join the worker threads
        for (int i = 0; i < num_threads; i++) {
                int res = pthread_join(worker_threads[i].thread_id, NULL);
                if (0 != res) {
                        perror("pthread_join failed");
                }
        }

        fprintf(stderr, "joined threads\n");

        // Clean up
        free(worker_threads);
        fprintf(stderr, "freed threads\n");
}

/* Allocate and initialize a Needleman-Wunsch alignment computation */
computation_t *
init_computation(char *s1, char *s2, int m, int k, int g)
{
        // Allocate for alignment computation instance
        fprintf(stderr, "allocating for computation\n");
        computation_t *C = (computation_t *)malloc(sizeof(computation_t));
        if (NULL == C) {
                perror("malloc failed");
                exit(1);
        }

        // We use an MxN table (M cols, N rows).  We add 1 to each of
        // the input strings' lengths to make room for the base
        // row/column (see init_table)
        fprintf(stderr, "beginning table dimension computation\n");
        int M = strlen(s1) + 1;
        fprintf(stderr, "strlen(s1)==%d\n", M);
        int N = strlen(s2) + 1;
        fprintf(stderr, "strlens computed\n");

        // Create and initialize the scores table
        C->scores_table = alloc_table(M, N);
        fprintf(stderr, "table allocated\n");
        init_table(C->scores_table, g, (num_threads > 1));
        fprintf(stderr, "table initialized\n");

        /* Alignment strings */
        C->top_string = s1;
        C->side_string = s2;

        /* Alignment scores/penalties */
        C->match_score = m;
        C->mismatch_penalty = k;
        C->gap_penalty = g;

        /* Total number of solutions founds */
        C->solution_count = 0;
        int res = pthread_rwlock_init(&C->solution_count_rwlock, NULL);
        if (0 != res) {
                perror("pthread_rwlock_init failed");
                exit(1);
        }

        return C;
}

void
free_computation(computation_t *C)
{
        free_table(C->scores_table, (num_threads > 1));
        int res = pthread_rwlock_destroy(&C->solution_count_rwlock);
        if (0 != res) {
                perror("pthread_rwlock_destroy failed");
                exit(1);
        }

        free(C);
}

void
needleman_wunsch(char *s1, char *s2, int m, int k, int g)
{
//        fprintf(stderr, "\nstarting needleman_wunsch()\n");
        // Initialize computation
        computation_t *C = init_computation(s1, s2, m, k, g);

        // Set global computation pointer
        comp = C;

        // Fill out table, i.e. compute the optimal score
//        fprintf(stderr, "beginning compute_table_scores()\n");
        compute_table_scores();

        // Walk the table.  Mark the optimal path if tflag is set, print
        // the results of the algorithm's run if sflag is set, and print
        // the aligned strings if qflag is not set
        if (qflag != 1 || sflag == 1 || tflag == 1) {
                mark_optimal_path_in_table(C);
        }

        // Print table if tflag was set
        if (tflag == 1) {
                // extra newline to separate the output sections
                if (qflag != 1 || sflag == 1) {
                        printf("\n");
                }
                print_table(C->scores_table, C->top_string, C->side_string, uflag);
        }

        // Clean up
        free_computation(C);
}

void
stdin_check_fgetc_err_and_eof(int eof_ok)
{
        /* Verify we didn't get an error */
        if (0 != ferror(stdin)) {
                perror("fgetc failed");
                exit(1);
        }

        /* Verify we're not already at the end of stdin */
        if (1 != eof_ok) {
                if (0 != feof(stdin)) {
                        fprintf(stderr, "error: file ended too early\n");
                        exit(1);
                }
        }
}

char *
read_sequence_from_stdin(int eof_ok)
{
        int c;
        int i = 0;
        int s_max = INPUT_STRING_BASE_SIZE;
        char *s = (char *)malloc(s_max * sizeof(char));
        while (EOF != (c = fgetc(stdin))) {
                if (isspace(c)) {
                        break;
                } else {
                        s[i] = (char)c;
                        i = i + 1;
                }

                /* If we're out of room, allocate more space */
                if (s_max == i) {
                        s = realloc(s, s_max + INPUT_STRING_BASE_SIZE);
                        if (NULL == s) {
                                perror("realloc failed");
                                exit(1);
                        } else {
                                s_max = s_max + INPUT_STRING_BASE_SIZE;
                        }
                }
        }

        stdin_check_fgetc_err_and_eof(eof_ok);

        s[i+1] = '\0';

        return s;
}

void
read_strings_from_stdin(char **s1, char **s2)
{
        /* Read the first string from stdin */
        char *X = read_sequence_from_stdin(0);

        /* Read out the rest of the whitespace */
        int c = (int)' ';
        while (isspace(c)) {
                c = fgetc(stdin);
        }

        stdin_check_fgetc_err_and_eof(0);

        /* Put the last character back, as it's part of the next sequence */
        int res = ungetc(c, stdin);
        if (EOF == res) {
                perror("ungetc failed");
                exit(1);
        }

        /* Read the second string from stdin */
        char *Y = read_sequence_from_stdin(1);

        /* Place X and Y in the global computation instance */
//        fprintf(stderr, "%s %s\n", X, Y);
        *s1 = X;
        *s2 = Y;
}

int
main(int argc, char **argv)
{
        // Strings to align
        char *s1;
        char *s2;

        // Scoring values
        int m, k, g;

        // Parse option flags
        extern char *optarg;
        extern int optind;
        int c;
        while ((c = getopt(argc, argv, "chlp:qstu")) != -1) {
                switch (c) {
                case 'c':
                        cflag = 1;
                        break;
                case 'h':
                        usage();
                        break;
                case 'l':
                        lflag = 1;
                        break;
                case 'p':
                        num_threads = atoi(optarg);
                        if (num_threads < 2) {
                                fprintf(stderr, "Error: num_threads == %d; " \
                                        "num_threads must be greater than 1\n",
                                        num_threads);
                        }
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

        // Parse operands
        if ((optind + 3) > argc) {
                usage();
        } else {
                read_strings_from_stdin(&s1, &s2);
                //fprintf(stderr, "%s %s\n", s1, s2);
                /* s1 = argv[optind + 0]; */
                /* s2 = argv[optind + 1]; */

                // Scoring values
                m = atoi(argv[optind + 0]);
                k = atoi(argv[optind + 1]);
                g = atoi(argv[optind + 2]);
        }

        // Solve
        needleman_wunsch(s1, s2, m, k, g);

        // Clean up
        free(s1);
        free(s2);

        return 0;
}
