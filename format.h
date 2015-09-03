// Copyright (c) 2015, Scott Cheloha.
// All rights reserved.

// See LICENSE for full copyright information.

// format.h -- ANSI color escape sequences for more colorful printing.
//             Prototypes for routines to set/unset formatting options.

#ifndef COLOR_H
#define COLOR_H

#define ANSI_CSI_OPEN   "\x1b["
#define ANSI_SGI_CLOSE  "m"
#define ANSI_BOLD       "1"
#define ANSI_FG_RED     "31"
#define ANSI_FG_GREEN   "32"
#define ANSI_FG_YELLOW  "33"
#define ANSI_FG_BLUE    "34"
#define ANSI_FG_MAGENTA "35"
#define ANSI_FG_CYAN    "36"
#define ANSI_FMT_RESET  "0"

#define BOLD_FMT ANSI_CSI_OPEN ANSI_BOLD ANSI_SGI_CLOSE
#define MATCH_FMT ANSI_CSI_OPEN ANSI_FG_BLUE /*";" ANSI_BOLD*/ ANSI_SGI_CLOSE
#define MISMATCH_FMT ANSI_CSI_OPEN ANSI_FG_RED /*";" ANSI_BOLD*/ ANSI_SGI_CLOSE
#define GAP_FMT ANSI_CSI_OPEN ANSI_FG_YELLOW /*";" ANSI_BOLD*/ ANSI_SGI_CLOSE
#define RESET_FMT ANSI_CSI_OPEN ANSI_FMT_RESET ANSI_SGI_CLOSE

// Flag for colorful output.  Externed in main file.
int cflag;

// Type describing the various formatting options we support
typedef enum {bold_fmt, match_fmt, mismatch_fmt, gap_fmt} fmt_t;

// Set the output formatting to any of the formats in fmt_t the enum
void set_fmt(fmt_t f);

// Reset the output fomatting
void reset_fmt();

#endif /* COLOR_H */
