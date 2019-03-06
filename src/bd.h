#ifndef GUARD_bd_h
#define GUARD_bd_h

#include "rng.h"
#include "info.h"
#include "tree.h"
#include "funs.h"

double bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
double bdhet(tree& x, xinfo& xi, dinfo& di, double* phi, pinfo& pi, RNG& gen);

#endif
