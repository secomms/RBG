#include <stdio.h>
#include "params.h"
#include "utils.h"
#include "time.h"
#include "rng.h"
#include "string.h"
#include "keccak.h"

void commit(const int sigma[code_length], const int non_sys_H[code_length-code_dimension][code_dimension], const int non_sys_base_G[M ][code_length - M], int v_list[N][code_length], int sigma_list[N][code_length], char commitments[N][2*SEEDBYTES], int accumulate_sigma[code_length], char c[HASHBYTES], char seedTree[2*N-1][SEEDBYTES]);

void second_commit(const int beta, const int e[code_length], const int sigma[code_length], const int non_sys_H[code_length-code_dimension][code_dimension], const int non_sys_base_G[M ][code_length - M], const int v_list[N][code_length], const int sigma_list[N][code_length], int accumulate_sigma[code_length], int tilde_e_list[N][code_length], char h[HASHBYTES]);

void print_v_and_sigma(int v_list[N][code_length], int sigma_list[N][code_length]);

void print_e_tilde(int tilde_e_list[N][code_length]);
