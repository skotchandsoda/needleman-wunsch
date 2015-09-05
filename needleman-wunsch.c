// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch - align two strings with the Needleman-Wunsch algorithm
//                    (see http://en.wikipedia.org/Needleman–Wunsch_algorithm)

#include <math.h>
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
int iflag = 0;
int qflag = 0;
int sflag = 0;
int tflag = 0;
int uflag = 0;

void
usage()
{
    fprintf(stderr,
            "usage: needleman-wunsch [-chiqstu] s1 s2 m k d\n"          \
            "Align two strings with the Needleman-Wunsch algorithm\n"   \
            "arguments:\n"                                              \
            "  s1   first string to align (the \"top string\")\n"       \
            "  s2   second string to align (the \"side string\")\n"     \
            "   m   match bonus\n"                                      \
            "   k   mismatch penalty\n"                                 \
            "   d   gap penalty\n"                                      \
            "options:\n"                                                \
            "  -c   color the output with ANSI escape sequences\n"      \
            "  -h   print this message\n"                               \
            "  -i   print match, mismatch, and gap counts for each "    \
            "alignment pair\n"                                          \
            "  -q   be quiet; don't print the optimal alignments "      \
            "(cancels the '-i' flag)\n"                                 \
            "  -s   print additional statistics about the "             \
            "alignment pairs\n"                                         \
            "  -t   print the scores table; only useful for "           \
            "shorter input strings\n"                                   \
            "  -u   use unicode arrows when printing the scores table\n" \
        );
    exit(1);
}

void
print_aligned_string_char(char *s1, char *s2, int n)
{
        // If the characters mismatch and they aren't gap characters,
        // mark the output if ANSI formatting is set (cflag == 1)
        if (s1[n] != s2[n] && s1[n] != GAP_CHAR && s2[n] != GAP_CHAR) {
                set_fmt(mismatch_fmt);
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
                        print_aligned_strings_and_counts(X, Y, n-1, iflag);
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
compute_table_scores(char *s1, char *s2, table_t *T, int m, int k, int d)
{
        int match = 0;
        int diag_val = 0;
        int gap_in_x = 0;
        int gap_in_y = 0;

        // If we're printing the table, initialize the largest value
        if (tflag == 1) {
                T->greatest_abs_val = 0;
        }

        // Mark each cell with a score and any relevant directional
        // information
        for (int i = 1; i < T->M; i++) {
                for (int j = 1; j < T->N; j++) {
                        // Compute the maximum score
                        diag_val = (s1[i-1] == s2[j-1] ? m : (-k));
                        match = T->cells[i-1][j-1].score + diag_val;
                        gap_in_x = T->cells[i][j-1].score - d;
                        gap_in_y = T->cells[i-1][j].score - d;
                        T->cells[i][j].score = max3(match, gap_in_x, gap_in_y);

                        // If we're printing the table, update the
                        // largest value
                        if (tflag == 1 &&
                            abs(T->cells[i][j].score) > T->greatest_abs_val) {
                                T->greatest_abs_val = abs(T->cells[i][j].score);
                        }

                        // Mark the relevant optimal paths.  Multiple
                        // optimal paths are possible in a single cell.
                        if (T->cells[i][j].score == match) {
                                T->cells[i][j].diag = 1;
                        }
                        if (T->cells[i][j].score == gap_in_x) {
                                T->cells[i][j].up = 1;
                        }
                        if (T->cells[i][j].score == gap_in_y) {
                                T->cells[i][j].left = 1;
                        }
                }
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
        while ((c = getopt(argc, argv, "chiqstu")) != -1) {
                switch (c) {
                case 'c':
                        cflag = 1;
                        break;
                case 'h':
                        usage();
                        break;
                case 'i':
                        iflag = 1;
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

        // Solve
        needleman_wunsch(s1, s2, m, k, d);

        return 0;
}
