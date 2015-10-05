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
 * format.h - Definitions of ANSI color escape sequences for more
 *            colorful printing.  Prototypes for functions implemented
 *            in format.c.
 */

#ifndef __FORMAT_H__
#define __FORMAT_H__

/* ASCII 27 is the ESC to start CSI formatting in the ANSI terminal
 * standard*/
#define ANSI_CSI_OPEN   "\x1b["

/* An 'm' ends SGI formatting. */
#define ANSI_SGI_CLOSE  "m"

/* Tricky: BOLD can mean bold and/or bright colors.  It varies from
 * emulator to emulator, so consistent coloring/font handling is
 * impossible without significant effort. */
#define ANSI_BOLD       "1"
#define ANSI_UNDERLINE  "4"
#define ANSI_FG_RED     "31"
#define ANSI_FG_GREEN   "32"
#define ANSI_FG_YELLOW  "33"
#define ANSI_FG_BLUE    "34"
#define ANSI_FG_MAGENTA "35"
#define ANSI_FG_CYAN    "36"
#define ANSI_FMT_RESET  "0"

/* The AIXTERM extensions are supported by some terminals.  If they were
 * supported by all, we could consistently have bright colors _without_
 * bold text, but alas: the world is cruel. */
#define AIX_FG_RED      "91"
#define AIX_FG_GREEN    "92"
#define AIX_FG_YELLOW   "93"
#define AIX_FG_BLUE     "34"
#define AIX_FG_MAGENTA  "95"
#define AIX_FG_CYAN     "96"

/* ANSI formatting we use when printing the table and aligned strings */
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
        ANSI_FG_YELLOW ";" ANSI_BOLD            \
        ANSI_SGI_CLOSE

#define MATCH_CHAR_FMT OPT_PATH_FMT

#define MISMATCH_CHAR_FMT MISMATCH_ARROW_FMT

#define GAP_CHAR_FMT GAP_ARROW_FMT

#define NO_OVERLAP_CHAR_FMT ""

#define RESET_FMT                               \
        ANSI_CSI_OPEN                           \
        ANSI_FMT_RESET                          \
        ANSI_SGI_CLOSE

/* Flag for ANSI-escape formatted output.  It is externed in
 * needleman-wunsch.c. */
int cflag;

/* fmt_t describes the various formatting options we support when
 * printing the aligned strings (the default behavior) and the table
 * (via the '-t' flag). */
typedef enum {
        top_string_fmt,
        side_string_fmt,
        opt_path_fmt,
        match_arrow_fmt,
        mismatch_arrow_fmt,
        gap_arrow_fmt,
        match_char_fmt,
        mismatch_char_fmt,
        no_overlap_char_fmt,
        gap_char_fmt
} fmt_t;

void set_fmt(fmt_t f);

void reset_fmt();

#endif /* __FORMAT_H__ */
