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
 * score-table.h - Definition for score table and score table-cell of
 *                 Needleman-Wunsch alignment computation.  Contains
 *                 prototypes for functions implemented in
 *                 score-table.c.
 */

#ifndef __TABLE_H__
#define __TABLE_H__

#include <pthread.h>

#include "walk-table.h"

/* arrow_t: A type describing directions in the scores table. */
/* typedef enum {left, up, diag} arrow_t; */

/* cell_t: A type describing a cell in the scores table. */
typedef struct score_table_cell {
        int score;
        /* int diag; */
        /* int left; */
        /* int up; */
        /* int diag_done; */
        /* int left_done; */
        /* int up_done; */
        /* arrow_t src_direction; */
        int match;
        /* int in_optimal_path; */
        int processed;
        pthread_mutex_t score_mutex;
        pthread_cond_t processed_cv;
} score_table_cell_t;

/* table_t: A type describing an MxN table of cells (i.e. matrix of
 *          cell_t). */
typedef struct score_table {
        int M;
        int N;
        score_table_cell_t **cells;
        int max_score;
        unsigned int max_abs_score;
        /* unsigned int branch_count; */
        /* pthread_rwlock_t branch_count_rwlock; */
} score_table_t;

/* Allocate an MxN table of score_table_cells */
score_table_t *alloc_score_table(int M, int N);

/* Print the score table */
void print_table(score_table_t *S,
                 walk_table_t *W,
                 char *s1,
                 char *s2,
                 int unicode);

/* Destroy a table of M score_table_cell_t pointers */
void free_score_table(score_table_t *S, unsigned int nthreads);

#endif /* __TABLE_H__ */
