#include "dkl.h"

#ifndef DKL_TIME_H
#define DKL_TIME_H

static struct timespec t;

double get_current_time() {
  clock_gettime(CLOCK_REALTIME, &t);
  return t.tv_sec + (t.tv_nsec / 1.0e9);
}

#endif /* DKL_TIME_H */
