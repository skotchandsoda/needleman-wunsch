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
#include "table.h"

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
                 int num_threads)
{
        /* We use an MxN table (M cols, N rows).  We add 1 to each of
           the input strings' lengths to make room for the base
           row and base column (see init_table()) */
        int M = strlen(s1) + 1;
        debug("Top string is %d characters long", M-1);
        int N = strlen(s2) + 1;
        debug("Side string is %d characters long", N-1);

        /* Create and initialize the scores table */
        debug("Allocating scores table");
        C->scores_table = alloc_table(M, N);
        debug("Initializing scores table");
        init_table(C->scores_table, d, (num_threads > 1));

        /* Alignment strings */
        C->top_string = s1;
        C->side_string = s2;

        /* Alignment scores/penalties */
        C->match_score = m;
        C->mismatch_penalty = k;
        C->indel_penalty = d;

        /* Total number of solutions found */
        C->solution_count = 0;
        int res = pthread_rwlock_init(&(C->solution_count_rwlock), NULL);
        check(0 == res, "pthread_rwlock_init failed");

        /* Number of threads to use in the scoring step */
        C->num_threads = num_threads;

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
        free_table(C->scores_table, (C->num_threads > 1));
        int res = pthread_rwlock_destroy(&C->solution_count_rwlock);
        check(0 == res, "pthread_rwlock_destroy failed");
        free(C);
}
