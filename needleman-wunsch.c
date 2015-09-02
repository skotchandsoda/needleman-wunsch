// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch -- align two strings with the Needleman-Wunsch algorithm

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "table.h"

int qflag = 0;
int tflag = 0;
int uflag = 0;

void
usage()
{
    fprintf(stderr,
            "usage: needleman-wunsch [-hqt] s1 s2 m k d\n"              \
            "Align two strings with the Needleman-Wunsch algorithm\n"   \
            "arguments:\n"                                              \
            "  s1   first string to align\n"                            \
            "  s2   second string to align\n"                           \
            "   m   score for matching characters\n"                    \
            "   k   penalty for mismatched characters\n"                \
            "   d   penalty for a gap\n"                                \
            "options:\n"                                                \
            "  -h   print this message\n"                               \
            "  -q   be quiet; don't print an optimal alignment for "    \
            "the input strings\n"                                       \
            "  -t   print the scores table "                            \
            "(pretty-print it with awk or column)\n"                    \
            "  -u   use unicode arrows when printing the  table\n"      \
        );
    exit(1);
}

int
print_aligned_strings(char *X, char *Y, int n, int print_count)
{
        if (print_count > 0) {
                printf("\n");
        }

        // Print the strings backwards
        for (int a = n-1; a > -1; a--) {
                printf("%c", X[a]);
        }
        printf("\n");
        for (int a = n-1; a > -1; a--) {
                printf("%c", Y[a]);
        }
        printf("\n");

        return print_count + 1;
}

int
reconstruct_optimal_solutions(char *s1,
                              char *s2,
                              table_t *T,
                              char *X,
                              char *Y,
                              int i,
                              int j,
                              int n,
                              int print_count)
{
        // Base case; we've reached the top-left corner of the table
        if ((T->cells[i][j].diag || T->cells[i][j].left || T->cells[i][j].up) == 0) {
                print_count = print_aligned_strings(X, Y, n, print_count);
        } else { // Recursive case; build optimal strings for each possible path
                if (T->cells[i][j].diag == 1) {
                        X[n] = s1[i-1];
                        Y[n] = s2[j-1];
                        print_count = reconstruct_optimal_solutions(s1, s2, T, X, Y, i-1, j-1, n+1, print_count);
                }
                if (T->cells[i][j].left == 1) {
                        X[n] = s1[i-1];
                        Y[n] = '-';
                        print_count = reconstruct_optimal_solutions(s1, s2, T, X, Y, i-1, j, n+1, print_count);
                }
                if (T->cells[i][j].up == 1) {
                        X[n] = '-';
                        Y[n] = s2[j-1];
                        print_count = reconstruct_optimal_solutions(s1, s2, T, X, Y, i, j-1, n+1, print_count);
                }
        }

        return print_count;
}


void
print_optimal_solutions(char *s1, char *s2, table_t *T)
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
        int i = T->M-1;  // position (x direction)
        int j = T->N-1;  // position (y direction)
        int n = 0;       // character count

        int print_count = reconstruct_optimal_solutions(s1, s2, T, X, Y, i, j, n, 0);
        printf("\n%d optimal alignment%s\n",
               print_count,
               (print_count > 1 ? "s" : ""));
}

int
max(int a, int b, int c)
{
        int m = a;
        (m < b) && (m = b); //these are not conditional statements.
        (m < c) && (m = c); //these are just boolean expressions.
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
                        T->cells[i][j].score = max(match, gap_in_x, gap_in_y);
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

        // Print table
        if (tflag == 1) {
                print_table(T, s1, s2, uflag);
        }

        // Print aligned strings
        if (qflag != 1) {
                print_optimal_solutions(s1, s2, T);
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
        while ((c = getopt(argc, argv, "hqtu")) != -1) {
                switch (c) {
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
