#ifndef DKL_H
#define DKL_H

#define _POSIX_C_SOURCE 200809L

//System
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#include "dkl_types.h"

#include "mkl_compat.h"

#include "dkl_macros.h"
#include "dkl_memory.h"
#include "dkl_const.h"
#include "dkl_error.h"
#include "dkl_sparse_matrix.h"
#include "dkl_time.h"
#include "dkl_args.h"
#include "dkl_wf.h"

extern int64_t dkl_errno;

#endif /* DKL_H */
