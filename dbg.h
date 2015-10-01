/*-
Copyright (c) 2010, Zed A. Shaw.  All rights reserved.
Copyright (c) 2015, Scott S. Cheloha.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the name of the Learn C The Hard Way, Zed A. Shaw, nor the
      names of its contributors may be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
-*/

#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

/* Global program name for logging */
char *prog;

/* Set the global program name */
void set_prog_name(char *name);


#ifdef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...)                                      \
        fprintf(stderr,                                    \
                "%s: debug: %s:%d: " M "\n",               \
                prog, __FILE__, __LINE__, ##__VA_ARGS__)
#endif

/* If there is no errno string, erase the colon printed to delineate the
 * errno string */
#define clean_errno() (errno == 0 ? "\b\b \b" : strerror(errno))

/* If NDEBUG is defined, we print error and warning messages without
 * __FILE__ and __LINE__ information.  Otherwise, debugging has been
 * enabled, so we print those details for the programmer. */

#ifdef NDEBUG
#define log_err(M, ...)                                                 \
        fprintf(stderr,                                                 \
                "%s: error: " M ": %s\n",                               \
                prog, ##__VA_ARGS__, clean_errno())
#else
#define log_err(M, ...)                                                 \
        fprintf(stderr,                                                 \
                "%s: error: %s:%d: " M ": %s\n",                        \
                prog, __FILE__, __LINE__, ##__VA_ARGS__, clean_errno())
#endif

#ifdef NDEBUG
#define log_warn(M, ...)                                                \
        fprintf(stderr,                                                 \
                "%s: warning:" M ": %s\n",                              \
                prog, ##__VA_ARGS__, clean_errno())
#else
#define log_warn(M, ...)                                                \
        fprintf(stderr,                                                 \
                "%s: warning: %s:%d: %s " M "\n",                       \
                prog, __FILE__, __LINE__, ##__VA_ARGS__, clean_errno())
#endif

#define log_info(M, ...)                                   \
        fprintf(stderr,                                    \
                "%s: " M "\n",                             \
                prog, ##__VA_ARGS__)

#define check(A, M, ...)                        \
        if (!(A)) {                             \
                log_err(M, ##__VA_ARGS__);      \
                exit(1);                        \
        }

#define sentinel(M, ...)                        \
                log_err(M, ##__VA_ARGS__);      \
                exit(1);

#define unreachable()                                           \
        sentinel("unreachable code executed; giving up")

#endif /* __dbg_h__ */
