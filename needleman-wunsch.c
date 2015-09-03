// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch -- align two strings with the Needleman-Wunsch algorithm

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "format.h"
#include "table.h"

#define GAP_CHAR '-'

// Global flags
// defined in format.h:
extern int cflag;

// defined here:
int qflag = 0;
int tflag = 0;
int uflag = 0;

void
usage()
{
    fprintf(stderr,
            "usage: needleman-wunsch [-chqtu] s1 s2 m k d\n"            \
            "Align two strings with the Needleman-Wunsch algorithm\n"   \
            "arguments:\n"                                              \
            "  s1   first string to align\n"                            \
            "  s2   second string to align\n"                           \
            "   m   score for matching characters\n"                    \
            "   k   penalty for mismatched characters\n"                \
            "   d   penalty for a gap\n"                                \
            "options:\n"                                                \
            "  -c   color the output with ANSI escape sequences\n"      \
            "  -h   print this message\n"                               \
            "  -q   be quiet; don't print an optimal alignment for "    \
            "the input strings\n"                                       \
            "  -t   print the scores table "                            \
            "(pretty-print it with awk or column)\n"                    \
            "  -u   use unicode arrows when printing the  table\n"      \
        );
    exit(1);
}

void
print_char_from_aligned_string(char *s1, char *s2, int n)
{
        // Establish the ANSI formatting if cflag is set
        if (s1[n] == s2[n]) {
                set_fmt(match_fmt);
        } else if (s1[n] == GAP_CHAR || s2[n] == GAP_CHAR) {
                set_fmt(gap_fmt);
        } else {
                set_fmt(mismatch_fmt);
        }

        // Print the character
        printf("%c", s1[n]);

        // Reset the formatting if cflag was set
        reset_fmt();
}

int
print_aligned_strings(char *X, char *Y, int n, int soln_count)
{
        if (soln_count > 0) {
                printf("\n");
        }

        // Print the strings backwards
        for (int i = n; i > -1; i--) {
                print_char_from_aligned_string(X, Y, i);
        }
        printf("\n");
        for (int i = n; i > -1; i--) {
                print_char_from_aligned_string(Y, X, i);
        }
        printf("\n");

        return soln_count + 1;
}

int
reconstruct_solutions(char *s1,
                      char *s2,
                      table_t *T,
                      char *X,
                      char *Y,
                      int i,
                      int j,
                      int n,
                      int soln_count)
{
        // Base case: We've reached the top-left corner of the table, i.e.
        //            the current cell_t has no direction bits set.
        if ((T->cells[i][j].diag || T->cells[i][j].left || T->cells[i][j].up) == 0) {
                if (qflag != 1) {
                        soln_count = print_aligned_strings(X, Y, n-1, soln_count);
                }
        }
        // Recursive case: We're not at the top-left corner of the table.
        //                 Make a recursive call for each direction bit set.
        else {
                // Mark the cell so that (later) the table printing
                // routine can mark the cell as part of an optimal path
                T->cells[i][j].in_optimal_path = 1;

                if (T->cells[i][j].diag == 1) {
                        X[n] = s1[i-1];
                        Y[n] = s2[j-1];
                        soln_count = reconstruct_solutions(s1, s2, T,
                                                            X, Y, i-1, j-1,
                                                            n+1, soln_count);
                }
                if (T->cells[i][j].left == 1) {
                        X[n] = s1[i-1];
                        Y[n] = '-';
                        soln_count = reconstruct_solutions(s1, s2, T,
                                                            X, Y, i-1, j,
                                                            n+1, soln_count);
                }
                if (T->cells[i][j].up == 1) {
                        X[n] = '-';
                        Y[n] = s2[j-1];
                        soln_count = reconstruct_solutions(s1, s2, T,
                                                            X, Y, i, j-1,
                                                            n+1, soln_count);
                }
        }

        return soln_count;
}


void
print_solutions(char *s1, char *s2, table_t *T)
{
        // Allocate buffers for the aligned strings.  In the worst case
        // they will be M+N characters long.
        int max_aligned_strlen = T->M + T->N;
        char *X = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        if (X == NULL) {
                perror("malloc failed");
                exit(1);
        }
        char *Y = (char *)malloc((max_aligned_strlen * sizeof(char)) + 1);
        if (Y == NULL) {
                perror("malloc failed");
                exit(1);
        }

        // We move through the table starting at the bottom-right corner
        int i = T->M - 1;  // position (x direction)
        int j = T->N - 1;  // position (y direction)
        int n = 0;         // character count

        int soln_count = reconstruct_solutions(s1, s2, T, X, Y, i, j, n, 0);
        if (qflag != 1) {
                printf("\n%d optimal alignment%s\n",
                       soln_count,
                       (soln_count > 1 ? "s" : ""));
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
compute_optimal_alignment(char *s1,
                          char *s2,
                          table_t *T,
                          int m,
                          int k,
                          int d)
{
        int match = 0;
        int diag_val = 0;
        int gap_in_x = 0;
        int gap_in_y = 0;
        T->greatest_abs_val = 0;
        for (int i = 1; i < T->M; i++) {
                for (int j = 1; j < T->N; j++) {
                        diag_val = (s1[i-1] == s2[j-1] ? m : (-k));
                        match = T->cells[i-1][j-1].score + diag_val;
                        gap_in_x = T->cells[i][j-1].score - d;
                        gap_in_y = T->cells[i-1][j].score - d;
                        T->cells[i][j].score = max3(match, gap_in_x, gap_in_y);
                        if (abs(T->cells[i][j].score) > T->greatest_abs_val) {
                                T->greatest_abs_val = abs(T->cells[i][j].score);
                        }
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

        // Fill out table, i.e. compute the optimal score and alignment
        compute_optimal_alignment(s1, s2, T, m, k, d);

        // Walk the table, mark the optimal path if tflag is set, and
        // print aligned strings if qflag is not set
        if (qflag != 1 || tflag == 1) {
                print_solutions(s1, s2, T);
        }
        // Print table
        if (tflag == 1) {
                // extra newline to separate the output sections
                if (qflag != 1) {
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

        // Parse options
        extern char *optarg;
        extern int optind;
        int c;
        while ((c = getopt(argc, argv, "chqtu")) != -1) {
                switch (c) {
                case 'c':
                        cflag = 1;
                        break;
                case 'h':
                        usage();
                        break;
                case 'q':
                        qflag = 1;
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

        // Parse positionals
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
