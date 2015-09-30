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
 * computation.h - Definition of computation instance and prototypes for
 *                 function implementations in computation.c.
 */

#ifndef __COMPUTATION_H__
#define __COMPUTATION_H__

#include "score-table.h"
#include "walk-table.h"

/* Instance of a Needleman-Wunsch alignment computation */
typedef struct computation {
        /* Sequences to align */
        char *top_string;
        char *side_string;

        /* Scoring parameters */
        int match_score;
        int mismatch_penalty;
        int indel_penalty;

        /* Score table */
        score_table_t *score_table;

        /* The walk_table maintains state during alignment
         * reconstruction */
        walk_table_t *walk_table;

        /* We track the solution count for summarization purposes as
         * well as for profiling the performance of our alignment
         * reconstruction routine(s). */
        unsigned int solution_count;
        pthread_rwlock_t solution_count_rwlock;

        /* Number of threads to execute in parallel when writing scores
         * to score_table (defined above). */
        unsigned int num_threads;

        /* A store for pointers to each of the worker threads we execute
         * in parallel when writing scores to score_table. */
        pthread_t *worker_threads;
} computation_t;

/*
 * Prototypes
 */

computation_t *alloc_computation();

void init_computation_tables(score_table_t *S,
                             walk_table_t *W,
                             int d,
                             unsigned int nthreads);

computation_t *init_computation(computation_t *C,
                                char *s1,
                                char *s2,
                                int m,
                                int k,
                                int d,
                                unsigned int nthreads);

void free_computation(computation_t *C);

void inc_solution_count(computation_t *C);

unsigned int get_solution_count(computation_t *C);

void print_summary(computation_t *C);

#endif /* __COMPUTATION_H__ */
