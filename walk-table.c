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
 * walk-table.c - The walk table is used to ... *drumroll* ... walk the
 *                score table (see score_table.c) in order to
 *                reconstruct alignments.  This file contains routines
 *                for allocating, initializing, and cleaning up the walk
 *                table in the Needleman Wunsch algorithm.
 */

#include <stdlib.h>

#include "dbg.h"
#include "walk-table.h"

walk_table_t *
alloc_walk_table(int M, int N)
{
        /* Allocate for the walk table */
        walk_table_t *W = (walk_table_t *)malloc(sizeof(walk_table_t));
        check(NULL != W, "malloc failed");

        W->M = M;
        W->N = N;

        /* Allocate top-level array for T->cells */
        W->cells = (walk_table_cell_t **)malloc(M * sizeof(walk_table_cell_t *));
        check(NULL != W->cells, "malloc failed");

        /* Allocate subarrays for T->cells */
        for (int i = 0; i < M; i++) {
                W->cells[i] = (walk_table_cell_t *)calloc(N, sizeof(walk_table_cell_t));
                check(NULL != W->cells[i], "malloc failed");
        }

        return W;
}

void
free_walk_table(walk_table_t *W, unsigned int nthreads)
{
        /* Destroy mutex and conditional variable objects */
        /* int res; */
        /* if (1 == multiple_threads) { */
        /*         for (int i = 0; i < W->M; i++) { */
        /*                 for (int j = 0; j < W->N; j++) { */
        /*                         res = pthread_mutex_destroy(&T->cells[i][j].score_mutex); */
        /*                         check(0 == res, "pthread_mutex_destroy failed"); */
        /*                         res = pthread_cond_destroy(&T->cells[i][j].processed_cv); */
        /*                         check(0 == res, "pthread_cond_destroy failed"); */
        /*                 } */
        /*         } */
        /* } */

        /* Free each subarray (column) of walk_table_cells */
        for (int i = 0; i < W->M; i++) {
                free(W->cells[i]);
        }

        /* Free the top-level array of walk_table_cell_t pointers */
        free(W->cells);

        int res = pthread_rwlock_destroy(&W->branch_count_rwlock);
        check(0 == res, "pthread_rwlock_destroy failed");

        /* Free the walk table itself */
        free(W);
}

void
inc_branch_count(walk_table_t *W, unsigned int nthreads)
{
        if (nthreads > 1) {
                pthread_rwlock_wrlock(&W->branch_count_rwlock);
        }

        W->branch_count = W->branch_count + 1;

        if (nthreads > 1) {
                pthread_rwlock_unlock(&W->branch_count_rwlock);
        }
}
