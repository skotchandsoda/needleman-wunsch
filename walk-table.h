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
 * walk-table.h - Definition of walk table and walk table cell of
 *                Needleman-Wunsch alignment reconstruction.  Contains
 *                prototypes for functions implemented in walk_table.c.
 */

#ifndef __WALK_TABLE_H__
#define __WALK_TABLE_H__

#include <pthread.h>

/* arrow_t: Directions in a walk_table_t. */
typedef enum {left, up, diag} arrow_t;

/* cell_t: Cell in a walk_table_t. */
typedef struct walk_table_cell {
        /* int score; */
        int diag;
        int left;
        int up;
        int diag_done;
        int left_done;
        int up_done;
        arrow_t src_direction;
        /* int match; */
        int in_optimal_path;
        /* int processed; */
        /* pthread_mutex_t score_mutex; */
        /* pthread_cond_t processed_cv; */
} walk_table_cell_t;

/* walk_table_t: An MxN table of walk_table_cells (i.e. matrix of
 *               walk_table_cell_t). */
typedef struct walk_table {
        int M;
        int N;
        walk_table_cell_t **cells;
        /* int greatest_abs_val; */
        unsigned int branch_count;
        pthread_rwlock_t branch_count_rwlock;
} walk_table_t;

/*
 * Prototypes
 */

walk_table_t *alloc_walk_table(int M, int N);

void free_walk_table(walk_table_t *W, unsigned int nthreads);

void inc_branch_count(walk_table_t *W, unsigned int nthreads);

#endif /* __WALK_TABLE_H__ */
