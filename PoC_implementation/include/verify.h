#ifndef VERIFY_H
#define VERIFY_H

#include <stdio.h>
#include <string.h>
#include "params.h"
#include "utils.h"
#include "rng.h"
#include "response.h"

int verify(const Resp *resp, const int i,
           const int non_sys_null_G[CODE_LENGTH_N - M][M],
           const int non_sys_base_G[M][CODE_LENGTH_N - M],
           const int non_sys_H[CODE_LENGTH_N -
                               CODE_DIMENSION_K][CODE_DIMENSION_K],
           const int beta, const int e[CODE_LENGTH_N],
           const int s[CODE_LENGTH_N - CODE_DIMENSION_K],
           const unsigned char c[HASHBYTES], const unsigned char h[HASHBYTES]);

#endif
