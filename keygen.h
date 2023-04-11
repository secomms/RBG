#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"
#include "utils.h"
#include "time.h"
#include "rng.h"
#include "string.h"

void key_gen(int non_sys_H[code_length-code_dimension][code_dimension], int non_sys_base_G[M ][code_length - M], int non_sys_null_G[code_length - M][M], int e[code_length], int sigma[code_length], int sigma_e[code_length], int s[code_length - code_dimension]);
void print_key_pair(int non_sys_H[code_length-code_dimension][code_dimension], int non_sys_base_G[M ][code_length - M], int non_sys_null_G[code_length - M][M], int e[code_length], int sigma[code_length], int sigma_e[code_length], int s[code_length - code_dimension]);
