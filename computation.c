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
 * computation.c - Allocation, initialization, and cleanup for an
 *                 instance of a Needleman-Wunsch alignment computation.
 */

#include "computation.h"
#include "dbg.h"
#include "stdlib.h"
#include "score-table.h"
#include "walk-table.h"

/*
 * alloc_computation()
 * -------------------
 * Allocate a Needleman-Wunsch alignment computation.
 */
computation_t *
alloc_computation()
{
        debug("Allocating for computation instance");
        computation_t *C = (computation_t *)malloc(sizeof(computation_t));
        check(NULL != C, "malloc failed");
        return C;
}

/*
 * init_computation_tables()
 * -------------------------
 * Initialize the scores table and the reference walk table for a
 * Needleman-Wunsch alignment computation.
 *
 * S: scores table to initialize
 * W: walk table to initialize
 * d: indel penalty (used to initialize the top-most row and left-most
 *    column with seed values for the scoring run
 * nthreads: number of threads we're using for this computation
 */
void
init_computation_tables(score_table_t *S, walk_table_t *W, int d, unsigned int nthreads)
{
        /* Initialize all mutex and condition variables. */
        int res;
        if (nthreads > 1) {
                for (int i = 0; i < S->M; i++) {
                        for (int j = 0; j < S->N; j++) {
                                res = pthread_mutex_init(&S->cells[i][j].score_mutex, NULL);
                                check(0 == res, "pthread_mutex_init failed");
                                res = pthread_cond_init (&S->cells[i][j].processed_cv, NULL);
                                check(0 == res, "pthread_cond_init failed");
                        }
                }
        }

        /* Initialize the largest value. */
        S->greatest_abs_val = 0;

        /* Initialize the table.  Cell (0,0) has a score of 0 and no
           optimal direction. */
        S->cells[0][0].score = 0;
        S->cells[0][0].processed = 1;
        W->cells[0][0].up_done = 1;
        W->cells[0][0].left_done = 1;
        W->cells[0][0].diag_done = 1;

        /* The rest of the topmost row has score i * (-d) and LEFT
         * direction. */
        for (int i = 1; i < S->M; i++) {
                S->cells[i][0].score = i * (-d);
                S->cells[i][0].processed = 1;
                W->cells[i][0].left = 1;
                W->cells[i][0].up_done = 1;
                W->cells[i][0].diag_done = 1;
        }

        /* The rest of the leftmost column has score j * (-d) and UP
         * direction. */
        for (int j = 1; j < S->N; j++) {
                S->cells[0][j].score = j * (-d);
                S->cells[0][j].processed = 1;
                W->cells[0][j].up = 1;
                W->cells[0][j].left_done = 1;
                W->cells[0][j].diag_done = 1;
        }

        W->branch_count = 0;
        res = pthread_rwlock_init(&(W->branch_count_rwlock), NULL);
        check(0 == res, "pthread_rwlock_init failed");
}

/* init_computation()
 * ------------------
 * Allocate and initialize a Needleman-Wunsch alignment computation.
 *
 * C:  computational instance to initialize
 * s1: top string, i.e. the string we are aligning against
 * s2: side string, i.e. the string we align against s1
 * m:  match bonus
 * k:  mismatch penalty
 * d:  indel, i.e. gap, penalty
 */
computation_t *
init_computation(computation_t *C,
                 char *s1,
                 char *s2,
                 int m,
                 int k,
                 int d,
                 unsigned int nthreads)
{
        /* We use an MxN table (M cols, N rows).  We add 1 to each of
           the input strings' lengths to make room for the base row and
           base column of scores (see init_score_table() in
           score-table.c). */
        int M = strlen(s1) + 1;
        debug("Top string is %d characters long", M-1);
        int N = strlen(s2) + 1;
        debug("Side string is %d characters long", N-1);

        /* Create and initialize the scores table */
        debug("Allocating score table");
        C->score_table = alloc_score_table(M, N);
        debug("Allocating walk table");
        C->walk_table = alloc_walk_table(M, N);
        debug("Initializing score and walk tables");
        init_computation_tables(C->score_table, C->walk_table, d, nthreads);

        /* Alignment strings */
        C->top_string = s1;
        C->side_string = s2;

        /* Alignment scores/penalties */
        C->match_score = m;
        C->mismatch_penalty = k;
        C->indel_penalty = d;

        /* Total number of solutions found */
        C->solution_count = 0;
        if (nthreads > 1) {
                int res = pthread_rwlock_init(&(C->solution_count_rwlock), NULL);
                check(0 == res, "pthread_rwlock_init failed");
        }

        /* Number of threads to use in the scoring step */
        C->num_threads = nthreads;

        return C;
}

/* free_computation()
 * ------------------
 * Clean up a Needleman-Wunsch alignment computation.
 *
 * C: the computation instance to clean up
 */
void
free_computation(computation_t *C)
{
        int res = 1;

        free_score_table(C->score_table, C->num_threads);
        free_walk_table(C->walk_table, C->num_threads);

        if (C->num_threads > 1) {
                res = pthread_rwlock_destroy(&C->solution_count_rwlock);
                check(0 == res, "pthread_rwlock_destroy failed");
        }

        free(C);
}

void
inc_solution_count(computation_t *C)
{
        if (C->num_threads > 1) {
                pthread_rwlock_wrlock(&C->solution_count_rwlock);
        }

        C->solution_count = C->solution_count + 1;

        if (C->num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }
}

unsigned int
get_solution_count(computation_t *C)
{
        if (C->num_threads > 1) {
                pthread_rwlock_rdlock(&C->solution_count_rwlock);
        }

        unsigned int count = C->solution_count;

        if (C->num_threads > 1) {
                pthread_rwlock_unlock(&C->solution_count_rwlock);
        }

        return count;
}

/*
 * print_summary()
 *
 *   Print details about the algorithm's run.  In particular, print the
 *   number of optimal alignments and the optimal alignment score.
 *
 *   C - computation instance to summarize
 */
void
print_summary(computation_t *C)
{
        unsigned int soln_count = get_solution_count(C);
        int max_col = C->score_table->M - 1;
        int max_row = C->score_table->N - 1;
        printf("%d optimal alignment%s\n",
               soln_count, (soln_count > 1 ? "s" : ""));
        printf("Optimal score is %-d\n",
               C->score_table->cells[max_col][max_row].score);
}
