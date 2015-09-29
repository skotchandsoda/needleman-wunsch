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
pthread_t *worker_threads;

/* Each worker thread has a thread id and a starting column.  The starting
   column is the first column in the scores table the thread will process.
   The thread then processes column start_col + num_threads, then
   start_col + 2*num_threads, etc. */
struct process_col_set_args {
        int start_col;    /* Column in scores_table to start at */
        computation_t *C; /* Needleman-Wunsch computation instance to */
                          /* process columns for */
};

#endif /* __NEEDLEMAN_WUNSCH_H__ */
