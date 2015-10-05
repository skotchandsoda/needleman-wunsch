# Copyright (c) 2015, Scott Cheloha.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#
#   3. Neither the name of the copyright holder nor the names of its
#      contributors may be used to endorse or promote products derived
#      from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PROGS = needleman-wunsch smith-waterman

NW_SRC = needleman-wunsch.c
NW_OBJ = needleman-wunsch.o
SW_SRC = smith-waterman.c
SW_OBJ = smith-waterman.o

SRC = score-table.c walk-table.c print-table.c \
      format.c dbg.c read-sequences.c computation.c
INC = $(SRC:.c=.h) needleman-wunsch.h
OBJ = ${SRC:.c=.o}

CFLAGS = -std=gnu99 -O3 -Wall -Wextra
LIB = -lpthread

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $<

all: CFLAGS += -DNDEBUG
all: $(PROGS)

debug: CFLAGS += -g
debug: $(PROGS)

needleman-wunsch: needleman-wunsch.o $(OBJ)
	$(CC) -o $@ needleman-wunsch.o $(OBJ) $(LIB)

smith-waterman: smith-waterman.o $(OBJ)
	$(CC) -o $@ smith-waterman.o $(OBJ) $(LIB)

# $(PROG): $(OBJ)
#	$(CC) -o $@ $(OBJ) $(LIB)

$(OBJ): $(INC)

.PHONY: clean
clean:
	rm -f $(OBJ) needleman-wunsch.o smith-waterman.o $(PROGS)
