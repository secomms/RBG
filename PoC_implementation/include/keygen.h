#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"
#include "utils.h"
#include "time.h"
#include "rng.h"
#include "string.h"

void key_gen(int non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K],
             int non_sys_base_G[M][CODE_LENGTH_N - M],
             int non_sys_null_G[CODE_LENGTH_N - M][M], int e[CODE_LENGTH_N],
             int sigma[CODE_LENGTH_N], int sigma_e[CODE_LENGTH_N],
             int s[CODE_LENGTH_N - CODE_DIMENSION_K]);
void print_key_pair(int
                    non_sys_H[CODE_LENGTH_N -
                              CODE_DIMENSION_K][CODE_DIMENSION_K],
                    int non_sys_base_G[M][CODE_LENGTH_N - M],
                    int non_sys_null_G[CODE_LENGTH_N - M][M],
                    int e[CODE_LENGTH_N], int sigma[CODE_LENGTH_N],
                    int sigma_e[CODE_LENGTH_N],
                    int s[CODE_LENGTH_N - CODE_DIMENSION_K]);
