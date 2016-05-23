#include <stdio.h>
#include "alloc_check.h"

// Define mallctl from jemalloc.
// - If we correctly link the jemalloc library, then this symbol will
//   be overridden.
// - If we do not link jemalloc, then, the printout below will be executed.
void mallctl(const char *name, void *oldp, size_t *oldlenp, void *newp, size_t newlen) {
    printf("Warning: jemalloc is not linked in.\n");
}
