// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// format.h -- ANSI color escape sequences for more colorful printing.
//             Prototypes for routines to set/unset formatting options.

#ifndef __FORMAT_H__
#define __FORMAT_H__

#define ANSI_CSI_OPEN   "\x1b[" // ASCII 27 is the ESC to start CSI formatting
#define ANSI_SGI_CLOSE  "m"     // 'm' ends SGI formatting
#define ANSI_BOLD       "1"
#define ANSI_UNDERLINE  "4"
#define ANSI_FG_RED     "31"
#define ANSI_FG_GREEN   "32"
#define ANSI_FG_YELLOW  "33"
#define ANSI_FG_BLUE    "34"
#define ANSI_FG_MAGENTA "35"
#define ANSI_FG_CYAN    "36"
#define ANSI_FMT_RESET  "0"

#define AIX_FG_RED      "91"
#define AIX_FG_GREEN    "92"
#define AIX_FG_YELLOW   "93"
#define AIX_FG_BLUE     "34"
#define AIX_FG_MAGENTA  "95"
#define AIX_FG_CYAN     "96"

// ANSI formatting we use when printing the table and aligned strings
#define TOP_STRING_FMT                          \
        ANSI_CSI_OPEN                           \
        ANSI_BOLD                               \
        ANSI_SGI_CLOSE

#define SIDE_STRING_FMT TOP_STRING_FMT

#define OPT_PATH_FMT                            \
        ANSI_CSI_OPEN                           \
        ANSI_FG_GREEN ";" ANSI_BOLD             \
        ANSI_SGI_CLOSE

#define MATCH_ARROW_FMT                         \
        ANSI_CSI_OPEN                           \
        ANSI_FG_CYAN ";" ANSI_BOLD              \
        ANSI_SGI_CLOSE

#define MISMATCH_ARROW_FMT                      \
        ANSI_CSI_OPEN                           \
        ANSI_FG_RED ";" ANSI_BOLD               \
        ANSI_SGI_CLOSE

#define GAP_ARROW_FMT                           \
        ANSI_CSI_OPEN                           \
        ANSI_FG_YELLOW ";" ANSI_BOLD           \
        ANSI_SGI_CLOSE

#define MATCH_CHAR_FMT ""

#define MISMATCH_CHAR_FMT MISMATCH_ARROW_FMT

#define GAP_CHAR_FMT ""

#define RESET_FMT                               \
        ANSI_CSI_OPEN                           \
        ANSI_FMT_RESET                          \
        ANSI_SGI_CLOSE

// Flag for ANSI-escape formatted output.  Externed in main file.
int cflag;

// Type describing the various formatting options we support when printing
// the aligned strings and the table
typedef enum {
        top_string_fmt,
        side_string_fmt,
        opt_path_fmt,
        match_arrow_fmt,
        mismatch_arrow_fmt,
        gap_arrow_fmt,
        match_char_fmt,
        mismatch_char_fmt,
        gap_char_fmt
} fmt_t;

// Set the output formatting to any of the formats described in the
// fmt_t enum
void set_fmt(fmt_t f);

// Reset the output fomatting
void reset_fmt();

#endif /* __FORMAT_H__ */
