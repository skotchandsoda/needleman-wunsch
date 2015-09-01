// Copyright (c) 2015 Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch -- align two strings with the Needleman-Wunsch algorithm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cell.h"
#include "table.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

int qflag = 0;
int tflag = 0;

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
        );
    exit(1);
}

void
print_aligned_strings(char *s1, char *s2, cell_t **table, int M, int N)
{
    char X[1024];
    char Y[1024];

    // Move through the table and modify the output strings as needed
    int i = M-1;
    int j = N-1;
    int n = 0;
    do {
        switch (table[i][j].ptr) {
        case DIAG:
            X[n] = s1[i-1];
            Y[n] = s2[j-1];
            i = i - 1;
            j = j - 1;
            break;
        case UP:
            X[n] = s1[i-1];
            Y[n] = '-';
            i = i - 1;
            break;
        case LEFT:
            X[n] = '-';
            Y[n] = s2[j-1];
            j = j - 1;
            break;
        case NONE:
            i = i - 1;
            j = j - 1;
            break;
        default:
            fprintf(stderr, "the impossible has happened; giving up\n");
            exit(1);
            break;
        }
        n = n + 1;
    } while (table[i][j].ptr != NONE);

    // Print the strings backwards
    for (int a = n-1; a > -1; a--) {
        printf("%c", X[a]);
    }
    printf("\n");
    for (int b = n-1; b > -1; b--) {
        printf("%c", Y[b]);
    }
    printf("\n");
}

void
compute_optimal_alignment(char *s1,
                          char *s2,
                          cell_t **table,
                          int M,
                          int N,
                          int m,
                          int k,
                          int d)
{
    int match = 0;
    int gap_in_x = 0;
    int gap_in_y = 0;
    for (int i = 1; i < M; i++) {
        for (int j = 1; j < N; j++) {
            match = table[i-1][j-1].score + (s1[i-1] == s2[j-1] ? m : (-k));
            gap_in_x = table[i-1][j].score - d;
            gap_in_y = table[i][j-1].score - d;
            table[i][j].score = MAX(MAX(match, gap_in_x),
                                    MAX(match, gap_in_y));
            if (table[i][j].score == match) {
                table[i][j].ptr = DIAG;
            } else if (table[i][j].score == gap_in_x) {
                table[i][j].ptr = UP;
            } else if (table[i][j].score == gap_in_y) {
                table[i][j].ptr = LEFT;
            } else {
                fprintf(stderr, "the impossible has happened; giving up\n");
                exit(1);
            }
        }
    }
}


void
needleman_wunsch(char *s1, char *s2, int m, int k, int d)
{
    // We use an MxN table (M cols, N rows).  We add 1 to each of the
    // input strings' lengths to make room for the base row/column (see
    // init_table)
    int M = strlen(s1) + 1;
    int N = strlen(s2) + 1;

    // Create and initialize the scores table
    cell_t **table = alloc_table(M, N);
    init_table(table, M, N, d);

    // Fill out table, i.e. compute the optimal score and alignment
    compute_optimal_alignment(s1, s2, table, M, N, m, k, d);

    // Print table
    if (tflag == 1) {
        print_table(table, M, N, s1, s2);
    }

    // Print aligned strings
    if (qflag != 1) {
        print_aligned_strings(s1, s2, table, M, N);
    }

    // Clean up
    free_table(table, M);
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
    while ((c = getopt(argc, argv, "hqt")) != -1) {
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
