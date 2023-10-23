#ifndef RESPONSE_H
#define RESPONSE_H

#include <string.h>
#include "params.h"

typedef struct {
   unsigned char c_i[HASHBYTES];
   int tilde_e_i[CODE_LENGTH_N];
   int sigma_1[CODE_LENGTH_N];
   unsigned char seedPath[LOGN][SEEDBYTES];
} Resp;

void pack_response(Resp *resp, const int i, const unsigned char *c_i,
                   const int *tilde_e_i, const int tilde_e_i_len,
                   const int *sigma, const int sigma_len,
                   const unsigned char seedPath[LOGN][SEEDBYTES]);

void unpack_response(const Resp *resp, const int i, unsigned char *c_i,
                     int *tilde_e_i, const int tilde_e_i_len,
                     int *sigma_1, const int sigma_1_len,
                     unsigned char seedPath[LOGN][SEEDBYTES]);

#endif
