#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include "string.h"
#include "params.h"

//generate random bytes
void sample_seed_and_salt(unsigned char stream[2 * SEEDBYTES]);

void sample_seed(unsigned char stream[SEEDBYTES]);

void sample_from_seed(const unsigned char seed[SEEDBYTES],
                      const int non_sys_base_G[M][CODE_LENGTH_N - M],
                      int sigma[CODE_LENGTH_N], int v[CODE_LENGTH_N],
                      unsigned char salt[SEEDBYTES]);

void sample_from_seed_only_v(const unsigned char seed[SEEDBYTES],
                             int v[CODE_LENGTH_N], unsigned char salt[SEEDBYTES]);

void gen_seed_tree(unsigned char seedTree[2 * N - 1][SEEDBYTES]);

void gen_seed_path(unsigned char seedPath[LOGN][SEEDBYTES],
                   const unsigned char seedTree[2 * N - 1][SEEDBYTES], const int i);

void genSeeds(unsigned char seeds[N][SEEDBYTES],
              const unsigned char seedPath[LOGN][SEEDBYTES], const int i);
