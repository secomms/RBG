#include "response.h"

/* This packing and unpacking of seeds is actually not very efficient */

void pack_response(Resp *resp, const int i, const unsigned char *c_i,
                   const int *tilde_e_i, const int tilde_e_i_len,
                   const int *sigma, const int sigma_len,
                   const unsigned char seedPath[LOGN][SEEDBYTES])
{

   memcpy(resp->c_i, c_i, HASHBYTES);
   memcpy(resp->tilde_e_i, tilde_e_i, tilde_e_i_len);

   if (i != 0) {
      memcpy(resp->sigma_1, sigma, sigma_len);
   }

   for (int k = 0; k < LOGN; k++) {
      memcpy(&resp->seedPath[k], seedPath[k], SEEDBYTES);
   }
}


void unpack_response(const Resp *resp, const int i, unsigned char *c_i,
                     int *tilde_e_i, const int tilde_e_i_len,
                     int *sigma_1, const int sigma_1_len,
                     unsigned char seedPath[LOGN][SEEDBYTES])
{
   memcpy(c_i, resp->c_i, HASHBYTES);
   memcpy(tilde_e_i, resp->tilde_e_i, tilde_e_i_len);

   if (i != 0) {
      memcpy(sigma_1, resp->sigma_1, sigma_1_len);
   }

   for (int k = 0; k < LOGN; k++) {
      memcpy(seedPath[k], &resp->seedPath[k], SEEDBYTES);
   }
}
