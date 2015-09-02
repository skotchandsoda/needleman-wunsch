# Copyright (C) 2015, Scott Cheloha
# All rights reserved.

# See LICENSE for full copyright information.

PROG = needleman-wunsch
SRC = needleman-wunsch.c table.c
INC = table.h
OBJ = ${SRC:.c=.o}
CFLAGS = -g

.SUFFIXES:
.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c $<

all: $(PROG)

$(PROG): $(OBJ)
	$(CC) -o $@ $(OBJ)

$(OBJ): $(INC)

clean:
	rm -f $(OBJ) $(PROG)
