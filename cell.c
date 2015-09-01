// Copyright (C) 2015 Scott Cheloha.
// All rights reserved.

// cell.c -- routine for printing a cell in a Needleman Wunsch table

#include <stdio.h>
#include <stdlib.h>

#include "cell.h"

void
print_cell(cell_t *c, int print_direction)
{
    printf("%d", c->score);

    if (print_direction) {
        printf("(");

        char dir;
        switch (c->ptr) {
        case DIAG:
            dir = '\\';
            break;
        case UP:
            dir = '^';
            break;
        case LEFT:
            dir = '<';
            break;
        case NONE:
            dir = 'X';
            break;
        default:
            fprintf(stderr, "the impossible has happened; giving up\n");
            exit(1);
            break;
        }

        printf("%c)", dir);
    }
}
