// Copyright (c) 2015, Scott Cheloha
// All rights reserved.

// See LICENSE for full copyright information.

// format.c -- Set/unset ANSI escape-sequence output formatting

#include <stdio.h>
#include <stdlib.h>

#include "format.h"

// Set the output formatting to any of the formats in fmt_t the enum
void
set_fmt(fmt_t f)
{
        if (cflag == 1) {
                switch (f) {
                case top_string_fmt:
                        printf(TOP_STRING_FMT);
                        break;
                case side_string_fmt:
                        printf(SIDE_STRING_FMT);
                        break;
                case opt_path_fmt:
                        printf(OPT_PATH_FMT);
                        break;
                case match_arrow_fmt:
                        printf(MATCH_ARROW_FMT);
                        break;
                case mismatch_arrow_fmt:
                        printf(MISMATCH_ARROW_FMT);
                        break;
                case gap_arrow_fmt:
                        printf(GAP_ARROW_FMT);
                        break;
                case match_char_fmt:
                        printf(MATCH_CHAR_FMT);
                        break;
                case mismatch_char_fmt:
                        printf(MISMATCH_CHAR_FMT);
                        break;
                case gap_char_fmt:
                        printf(GAP_CHAR_FMT);
                        break;
                default:
                        fprintf(stderr, "the impossible has happened; giving up\n");
                        exit(1);
                        break;
                }
        }
}

// Reset the output fomatting
void reset_fmt()
{
        if (cflag == 1) {
                printf(RESET_FMT);
        }
}
