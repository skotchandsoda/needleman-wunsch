// Copyright (c) 2015, Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// format.h -- ANSI color escape sequences for more colorful printing.
//             Prototypes for routines to set/unset formatting options.

#ifndef COLOR_H
#define COLOR_H

/* #define ANSI_CSI_OPEN   "\x1b[" */
/* #define ANSI_SGI_CLOSE  "m" */
/* #define ANSI_BOLD       "1" */
/* #define ANSI_FG_RED     "31" */
/* #define ANSI_FG_GREEN   "32" */
/* #define ANSI_FG_YELLOW  "33" */
/* #define ANSI_FG_BLUE    "34" */
/* #define ANSI_FG_MAGENTA "35" */
/* #define ANSI_FG_CYAN    "36" */
/* #define ANSI_FMT_RESET  "0" */

#define SCORE_FMT ""
#define MATCH_FMT "\x1b[34;1m"
#define MISMATCH_FMT "\x1b[31;1m"
#define GAP_FMT "\x1b[30;1m"
#define RESET_FMT "\x1b[0m"

// Flag for colorful output.  Externed in main file.
int cflag;

// Type describing the various formatting options we support
typedef enum {score_fmt, match_fmt, mismatch_fmt, gap_fmt} fmt_t;

// Set the output formatting to any of the formats in fmt_t the enum
void set_fmt(fmt_t f);

// Reset the output fomatting
void reset_fmt();

#endif /* COLOR_H */
