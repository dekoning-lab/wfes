#include "dkl.h"

#ifndef DKL_TIME_H
#define DKL_TIME_H

static struct timespec t;

double get_current_time() {
#ifdef __MACH__
  clock_serv_t cclock;
  mach_timespec_t mt;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mt);
  mach_port_deallocate(mach_task_self(), cclock);
  t.tv_sec = mt.tv_sec;
  t.tv_nsec = mt.tv_nsec;
#else
  clock_gettime(CLOCK_REALTIME, &t);
#endif
  return t.tv_sec + (t.tv_nsec / 1.0e9);
}

#endif /* DKL_TIME_H */
