// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch - align two strings with the Needleman-Wunsch algorithm
//                    (see http://en.wikipedia.org/Needleman–Wunsch_algorithm)

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "format.h"
#include "table.h"

#define GAP_CHAR '-'

// Global flags defined in format.h
extern int cflag;

// Global flags defined in needleman-wunsch.c (this file)
int lflag = 0;
int qflag = 0;
int sflag = 0;
int tflag = 0;
int uflag = 0;

int num_threads = 1;

void
usage()
{
        fprintf(stderr,"\
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
  -l   print match, mismatch, and gap counts for each alignment pair\n\
  -p num_threads\n\
       parallelize the computation with 'num_threads' threads (must be >1)\n\
  -q   be quiet and don't print the alignmened strings (cancels the '-i' flag)\n\
  -s   summarize the algorithm's run\n\
  -t   print the scores table; probably only useful for shorter input strings\n\
  -u   use unicode arrows when printing the scores table\n"
);
        exit(1);
}

void
print_aligned_string_char(char *s1, char *s2, int n)
{
        // Format the output character as defined in format.h
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

        // Print the character
        printf("%c", s1[n]);

        reset_fmt();
}

void
print_aligned_strings_and_counts(char *X, char *Y, int n, int print_counts)
{
        int match_count = 0;
        int mismatch_count = 0;
        int gap_count = 0;

        // Print the strings backwards
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

int
walk_table_recursively(char *s1,
                       char *s2,
                       table_t *T,
                       char *X,
                       char *Y,
                       int i,
                       int j,
                       int n,
                       int soln_count)
{
        // Mark the cell so that (later) the table printing
        // routine can mark the cell as part of an optimal path
        if (tflag == 1) {
                T->cells[i][j].in_optimal_path = 1;
        }

        // Base case: We've reached the top-left corner of the table, i.e.
        //            the current cell_t has no direction bits set.
        if ((T->cells[i][j].diag || T->cells[i][j].left || T->cells[i][j].up) == 0) {
                if (qflag != 1) {
                        if (soln_count > 0) {
                                printf("\n");
                        }
                        print_aligned_strings_and_counts(X, Y, n-1, lflag);
                }
                soln_count = soln_count + 1;
        }
        // Recursive case: We're not at the top-left corner of the table.
        //                 Make a recursive call for each direction bit set.
        else {
                if (T->cells[i][j].diag == 1) {
                        if (qflag != 1) {
                                X[n] = s1[i-1];
                                Y[n] = s2[j-1];
                        }
                        soln_count = walk_table_recursively(s1, s2, T,
                                                            X, Y, i-1, j-1,
                                                            n+1, soln_count);
                }
                if (T->cells[i][j].left == 1) {
                        if (qflag != 1) {
                                X[n] = s1[i-1];
                                Y[n] = '-';
                        }
                        soln_count = walk_table_recursively(s1, s2, T,
                                                            X, Y, i-1, j,
                                                            n+1, soln_count);
                }
                if (T->cells[i][j].up == 1) {
                        if (qflag != 1) {
                                X[n] = '-';
                                Y[n] = s2[j-1];
                        }
                        soln_count = walk_table_recursively(s1, s2, T,
                                                            X, Y, i, j-1,
                                                            n+1, soln_count);
                }
        }

        return soln_count;
}


void
mark_optimal_path_in_table(char *s1, char *s2, table_t *T)
{
        int max_aligned_strlen;
        char *X;
        char *Y;

        // If the qflag wasn't set, allocate buffers for printing the
        // optimally aligned strings.  In the worst case they will be
        // M+N characters long.
        if (qflag != 1) {
                max_aligned_strlen = T->M + T->N;
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

        // We move through the table starting at the bottom-right corner
        int i = T->M - 1;  // position (x direction)
        int j = T->N - 1;  // position (y direction)
        int n = 0;         // character count

        // Walk the table starting at the bottom-right corner recursively,
        // marking cells in the optimal path and counting the total possible
        // optimal solutions (alignments)
        int soln_count = walk_table_recursively(s1, s2, T, X, Y, i, j, n, 0);

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
                printf("%d optimal alignment%s\n",
                       soln_count, (soln_count > 1 ? "s" : ""));
                printf("Optimal score is %-d\n", T->cells[i][j].score);
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
process_cell(table_t *T, int col, int row, char *s1, char *s2, int m, int k, int d)
{
        // Cell we want to compute the score for
        cell_t *target_cell = &T->cells[col][row];

        // Cells we'll use to compute that score
        cell_t *up_cell   = &T->cells[col][row-1];
        cell_t *left_cell = &T->cells[col-1][row];
        cell_t *diag_cell = &T->cells[col-1][row-1];

        // Candidate scores
        int up_score = up_cell->score - d;
        int left_score = left_cell->score - d;
        int diag_score = 0;
        if (s1[col-1] == s2[row-1]) {
                diag_score = diag_cell->score + m;
                target_cell->match = 1;
        } else {
                diag_score = diag_cell->score - k;
                target_cell->match = 0;
        }

        // The current cell's score is the max of the three candidate scores
        target_cell->score = max3(up_score, left_score, diag_score);

        // Mark the relevant optimal paths.  Provided that a path's
        // score is equal to the target cell's score, i.e. the maximum
        // of the three candidate scores, it is an optimal path.
        if (target_cell->score == diag_score) {
                target_cell->diag = 1;
        }
        if (target_cell->score == up_score) {
                target_cell->up = 1;
        }
        if (target_cell->score == left_score) {
                target_cell->left = 1;
        }

        target_cell->processed = 1;
}

void
process_column(table_t *T, int col, char *s1, char *s2, int m, int k, int d)
{
        // Compute the score for each cell in the column
        for (int row = 1; row < T->N; row++) {
                // Compute the cell's score
                process_cell(T, col, row, s1, s2, m, k, d);

                // If we're printing the table and the absolute value of
                // the current cell's score is greater than the one
                // marked in the table, update the largest value
                int current_abs_score = abs(T->cells[col][row].score);
                if (tflag == 1 && current_abs_score > T->greatest_abs_val) {
                        T->greatest_abs_val = current_abs_score;
                }
        }
}

void
compute_table_scores(char *s1, char *s2, table_t *T, int m, int k, int d)
{
        // If we're printing the table, initialize the largest value
        if (tflag == 1) {
                T->greatest_abs_val = 0;
        }

        // Mark each cell with a score and any relevant directional
        // information
        for (int col = 1; col < T->M; col++) {
                process_column(T, col, s1, s2, m, k, d);
        }
}


void
needleman_wunsch(char *s1, char *s2, int m, int k, int d)
{
        // We use an MxN table (M cols, N rows).  We add 1 to each of
        // the input strings' lengths to make room for the base
        // row/column (see init_table)
        int M = strlen(s1) + 1;
        int N = strlen(s2) + 1;

        // Create and initialize the scores table
        table_t *T = alloc_table(M, N);
        init_table(T, d);

        // Fill out table, i.e. compute the optimal score
        compute_table_scores(s1, s2, T, m, k, d);

        // Walk the table.  Mark the optimal path if tflag is set, print
        // the results of the algorithm's run if sflag is set, and print
        // the aligned strings if qflag is not set
        if (qflag != 1 || sflag == 1 || tflag == 1) {
                mark_optimal_path_in_table(s1, s2, T);
        }

        // Print table if tflag was set
        if (tflag == 1) {
                // extra newline to separate the output sections
                if (qflag != 1 || sflag == 1) {
                        printf("\n");
                }
                print_table(T, s1, s2, uflag);
        }

        // Clean up
        free_table(T);
}

int
main(int argc, char **argv)
{
        // Strings to align
        char *s1;
        char *s2;

        // Scoring values
        int m, k, d;

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
                                usage();
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

        // Parse positional arguments
        if ((optind + 5) > argc) {
                usage();
        } else {
                s1 = argv[optind + 0];
                s2 = argv[optind + 1];

                // Scoring values
                m = atoi(argv[optind + 2]);
                k = atoi(argv[optind + 3]);
                d = atoi(argv[optind + 4]);
        }

        fprintf(stderr, "num_threads == %d\n", num_threads);

        // Solve
        needleman_wunsch(s1, s2, m, k, d);

        return 0;
}
