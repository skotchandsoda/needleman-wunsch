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
 * table.h - Definition for scores table and scores table cell of
 *           Needleman-Wunsch alignment computation.  Also contains
 *           prototypes for functions implemented in table.c.
 */

#ifndef __TABLE_H__
#define __TABLE_H__

#include <pthread.h>

/* arrow_t: A type describing directions in the scores table. */
typedef enum {left, up, diag} arrow_t;

/* cell_t: A type describing a cell in the scores table. */
typedef struct cell {
        int score;
        int diag;
        int left;
        int up;
        int diag_done;
        int left_done;
        int up_done;
        arrow_t src_direction;
        int match;
        int in_optimal_path;
        int processed;
        pthread_mutex_t score_mutex;
        pthread_cond_t processed_cv;
} cell_t;

/* table_t: A type describing an MxN table of cells (i.e. matrix of
 *          cell_t). */
typedef struct table {
        int M;
        int N;
        cell_t **cells;
        int greatest_abs_val;
        unsigned int branch_count;
        pthread_rwlock_t branch_count_rwlock;
} table_t;

/* Allocate a MxN table of cells */
table_t *alloc_table(int M, int N);

/* Destroy a table of M cell_t pointers */
void free_table(table_t *T, int multiple_threads);

/* Initialize the score table */
void init_table(table_t *T, int d, int multiple_threads);

/* Print the score table */
void print_table(table_t *T, char *s1, char *s2, int unicode);

#endif /* __TABLE_H__ */
