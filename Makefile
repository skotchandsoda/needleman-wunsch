# Copyright (C) 2015 Scott Cheloha.
# All rights reserved.

PROG := needleman-wunsch
SRC := needleman-wunsch.c table.c cell.c
OPTS :=

all: $(PROG)

debug: OPTS += -g
debug $(PROG): $(SRC)
	$(CC) $(OPTS) -o $(PROG) $(SRC)

clean:
	rm -rf *.o *.dSYM $(PROG)
