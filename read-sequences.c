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

/* read-sequences.c - Read sequences from a stream into memory. */

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "dbg.h"
#include "read-sequences.h"

static char *
read_sequence_from_stream(FILE *in)
{
        int c;
        int i = 0;
        int seq_max = INPUT_STRING_BUF_SIZE;
        char *seq = (char *)malloc(seq_max * sizeof(char));
        check(NULL != seq, "malloc failed");

        /* Read characters from stdin into seq until we hit whitespace */
        while (EOF != (c = fgetc(in)) && !isspace(c)) {
                seq[i] = (char)c;
                i = i + 1;

                /* If we're out of room, try to allocate more space */
                if (seq_max == i) {
                        seq = realloc(seq, seq_max + INPUT_STRING_BUF_SIZE);
                        check(NULL != seq, "realloc failed");
                        seq_max = seq_max + INPUT_STRING_BUF_SIZE;
                }
        }

        /* Null-terminate the buffer by hand.  Note that i is the last
         * index we never wrote a character to. */
        seq[i] = '\0';

        return seq;
}

/*
 * discard_whitespace_in_stream()
 *
 *   Read in characters until isspace(3) returns false.
 *
 *   in - stream to read from
 *
 *   return - first non-whitespace character read from stream
 */
static int
discard_whitespace_in_stream(FILE *in)
{
        int c = (int)' ';
        while (isspace(c))
                c = fgetc(in);

        return c;
}

void
read_two_sequences_from_stream(char **s1, char **s2, FILE *in)
{
        /* Read the first string from the input stream */
        char *X = read_sequence_from_stream(in);

        /* Verify we didn't get an error */
        check(0 == ferror(in), "fgetc failed while reading input string 1");

        /* Verify we didn't already read EOF */
        check(0 == feof(in), "got EOF too early while reading input strings; " \
              "expected at least two input strings but found only 1");

        /* Discard the whitespace between the two input strings. */
        int c = discard_whitespace_in_stream(in);

        /* Verify we didn't get an error */
        check(0 == ferror(in), "fgetc failed while clearing whitespace " \
              "between input strings");

        /* Verify we didn't already read EOF */
        check(0 == feof(in), "got EOF too early when reading input strings; " \
              "expected at least two input strings but found only 1");

        /* Put the last character back, as it's part of the next sequence */
        int res = ungetc(c, in);
        check(EOF != res, "ungetc failed");

        /* Read the second string from the input stream */
        char *Y = read_sequence_from_stream(in);

        /* Check that we didn't get an error */
        check(0 == ferror(in), "fgetc failed while reading input string 2");

        /* Make X & Y available to the caller */
        *s1 = X;
        *s2 = Y;
}
