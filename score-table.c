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
 * score-table.c - Routines for allocating, initializing, and printing
 *                 the scoring table in the Needleman Wunsch algorithm.
 */

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "dbg.h"
#include "format.h"
#include "score-table.h"
#include "walk-table.h"


/*
 * alloc_score_table()
 *
 *   Allocate a score_table_t instance.
 *
 *   M - number of columns in the table
 *
 *   N - number of rows in the table
 *
 *   return - allocated score_table_t pointer with an allocated MxN
 *            score_table_cell_t matrix
 */
score_table_t *
alloc_score_table(int M, int N)
{
        /* Allocate for the scores table */
        score_table_t *S = (score_table_t *)malloc(sizeof(score_table_t));
        check(NULL != S, "malloc failed");

        S->M = M;
        S->N = N;

        /* Allocate top-level array for S->cells */
        S->cells = (score_table_cell_t **)malloc(M * sizeof(score_table_cell_t *));
        check(NULL != S->cells, "malloc failed");

        /* Allocate subarrays for S->cells */
        for (int i = 0; i < M; i++) {
                S->cells[i] = (score_table_cell_t *)calloc(N, sizeof(score_table_cell_t));
                check(NULL != S->cells[i], "calloc failed");
        }

        return S;
}

void
free_score_table(score_table_t *S, unsigned int nthreads)
{
        /* Destroy mutex and conditional variable objects */
        int res;
        if (nthreads > 1) {
                for (int i = 0; i < S->M; i++) {
                        for (int j = 0; j < S->N; j++) {
                                res = pthread_mutex_destroy(
                                        &S->cells[i][j].score_mutex);
                                check(0 == res, "pthread_mutex_destroy failed");
                                res = pthread_cond_destroy(
                                        &S->cells[i][j].processed_cv);
                                check(0 == res, "pthread_cond_destroy failed");
                        }
                }
        }

        /* Free each subarray (column) of cells */
        for (int i = 0; i < S->M; i++) {
                free(S->cells[i]);
        }

        /* Free the top-level array of cell pointers */
        free(S->cells);

        /* Free the scores table itself */
        free(S);
}
