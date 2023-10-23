#include <stdio.h>
#include "params.h"
#include "utils.h"
#include "time.h"
#include "rng.h"
#include "string.h"

void commit(const int sigma[CODE_LENGTH_N],
        const int non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K],
        const int non_sys_base_G[M][CODE_LENGTH_N - M],
        int v_list[N][CODE_LENGTH_N], int sigma_list[N][CODE_LENGTH_N],
        unsigned char commitments[N][2 * SEEDBYTES],
        int accumulate_sigma[CODE_LENGTH_N], 
        unsigned char c[HASHBYTES],
        unsigned char seedTree[2 * N - 1][SEEDBYTES]);

void second_commit(const int beta,
                   const int e[CODE_LENGTH_N],
                   const int sigma[CODE_LENGTH_N],
                   const int non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K],
                   const int non_sys_base_G[M][CODE_LENGTH_N - M],
                   const int v_list[N][CODE_LENGTH_N],
                   const int sigma_list[N][CODE_LENGTH_N],
                   int accumulate_sigma[CODE_LENGTH_N],
                   int tilde_e_list[N][CODE_LENGTH_N],
                   unsigned char h[HASHBYTES]);

void print_v_and_sigma(int v_list[N][CODE_LENGTH_N],
                       int sigma_list[N][CODE_LENGTH_N]);

void print_e_tilde(int tilde_e_list[N][CODE_LENGTH_N]);
