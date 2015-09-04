// Copyright (c) 2015, Scott Cheloha.
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
                case score_fmt:
                        printf(SCORE_FMT);
                        break;
                case match_fmt:
                        printf(MATCH_FMT);
                        break;
                case mismatch_fmt:
                        printf(MISMATCH_FMT);
                        break;
                case gap_fmt:
                        printf(GAP_FMT);
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
