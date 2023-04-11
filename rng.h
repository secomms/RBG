#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include "string.h"
#include "params.h"
#include "keccak.h"

//generate random bytes
void sample_seed_and_salt (char stream[2*SEEDBYTES]);

void sample_seed (char stream[SEEDBYTES]);

void sample_from_seed(const char seed[SEEDBYTES], const int non_sys_base_G[M][code_length - M], int sigma[code_length], int v[code_length],  char salt[SEEDBYTES]);

void sample_from_seed_only_v(const char seed[SEEDBYTES], int v[code_length],  char salt[SEEDBYTES]);

void gen_seed_tree(char seedTree[2*N-1][SEEDBYTES]);

void gen_seed_path(char seedPath[LOGN][SEEDBYTES], const char seedTree[2*N-1][SEEDBYTES], const int i);

void genSeeds(char seeds[N][SEEDBYTES], const char seedPath[LOGN][SEEDBYTES], const int i);
