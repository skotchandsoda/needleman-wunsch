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

#include "table.h"

/* Instance of a Needleman-Wunsch alignment computation */
typedef struct computation {
        char *top_string;
        char *side_string;
        int match_score;
        int mismatch_penalty;
        int indel_penalty;
        table_t *scores_table;
        unsigned int solution_count;
        pthread_rwlock_t solution_count_rwlock;
        int num_threads;
        pthread_t *worker_threads;
} computation_t;

/*
 * Prototypes
 */

computation_t *alloc_computation();

computation_t *init_computation(computation_t *C,
                                char *s1,
                                char *s2,
                                int m,
                                int k,
                                int d,
                                int num_threads);

void free_computation(computation_t *C);

#endif /* __COMPUTATION_H__ */
