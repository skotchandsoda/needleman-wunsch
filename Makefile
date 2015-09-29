# Copyright (c) 2015, Scott Cheloha.
# All rights reserved.

# See LICENSE for full copyright information.

PROG = needleman-wunsch
SRC = needleman-wunsch.c table.c format.c dbg.c read-sequences.c computation.c
INC = needleman-wunsch.h table.h format.h dbg.h read-sequences.h computation.h
OBJ = ${SRC:.c=.o}
CFLAGS = -g -O2
LIB = -lpthread

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $<

all: $(PROG)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

$(OBJ): $(INC)

clean:
	rm -f $(OBJ) $(PROG)
