// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// needleman-wunsch.h - prototypes and global flags for needleman-wunsch.c

#ifndef __NEEDLEMAN_WUNSCH_H__
#define __NEEDLEMAN_WUNSCH_H__

#include <pthread.h>

#include "table.h"

/*
 * Global flags affecting program logic
 */

int lflag = 0;
int qflag = 0;
int sflag = 0;
int tflag = 0;
int uflag = 0;

/*
 * Threading globals
 */

int num_threads = 1;

/* Each worker thread has a thread id and a starting column.  The starting
   column is the first column in the scores table the thread will process.
   The thread then processes column start_col + num_threads, then
   start_col + 2*num_threads, etc. */
typedef struct worker_thread {
        pthread_t thread_id;
        int start_col; // Column to process first in the scores table
} worker_thread_t;

/* Global array of struct worker_thread */
worker_thread_t *worker_threads;

/*
 * Needleman-Wunsch computation instance definitions/globals
 */

/* Instance of a Needleman-Wunsch alignment computation */
typedef struct computation {
        char *top_string;
        char *side_string;
        int match_score;
        int mismatch_penalty;
        int gap_penalty;
        table_t *scores_table;
        unsigned int solution_count;
        pthread_rwlock_t solution_count_rwlock;
} computation_t;

/* Global computation instance for this process */
computation_t *comp;

#endif /* __NEEDLEMAN_WUNSCH_H__ */
