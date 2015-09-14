#include "dbg.h"

/* Set the program's name for logging.
   Note:
     1. 'prog' must be declared in dbg.h.
     2. 'name' needs to remain allocated throughout the program's run. */
void
set_prog_name(char *name)
{
        /* Increment the pointer in the program is being run outside of $PATH */
        if (name[0] == '.' && name[1] == '/') {
                name = name + 2;
        }
        prog = name;
}
