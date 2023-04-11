#ifndef VERIFY_H
#define VERIFY_H

#include <stdio.h>
#include <string.h>
#include "params.h"
#include "utils.h"
#include "rng.h"
#include "response.h"

int verify( const Resp *resp, const int i, const int non_sys_null_G[code_length - M][M], const int non_sys_base_G[M][code_length - M], 
            const int non_sys_H[code_length-code_dimension][code_dimension], const int beta, const int e[code_length], const int s[code_length - code_dimension],
            const char c[HASHBYTES], const char h[HASHBYTES]);

#endif
