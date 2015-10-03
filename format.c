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

/* format.c - Set or unset ANSI escape-sequence output formatting. */

#include <stdio.h>
#include <stdlib.h>

#include "dbg.h"
#include "format.h"

/*
 * set_fmt()
 *
 *   Set the output formatting to any of the formats in the fmt_t enum.
 *   See format.h for the #defines used here.
 */
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
                        unreachable();
                        break;
                }
        }
}

/*
 * reset_fmt()
 *
 *   Reset the output fomatting (i.e. reset colors/bolding to normal) on
 *   the standard output.
 */
void reset_fmt()
{
        if (cflag == 1) {
                printf(RESET_FMT);
        }
}
