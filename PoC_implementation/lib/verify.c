#include "verify.h"
#include "sha3.h"
int verify(const Resp *resp, const int i,
           const int non_sys_null_G[CODE_LENGTH_N - M][M],
           const int non_sys_base_G[M][CODE_LENGTH_N - M],
           const int non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K],
           const int beta, const int e[CODE_LENGTH_N],
           const int s[CODE_LENGTH_N - CODE_DIMENSION_K],
           const unsigned char c[HASHBYTES], const unsigned char h[HASHBYTES])
{
   unsigned char commitments[N][HASHBYTES];
   int tilde_e_list[N][CODE_LENGTH_N];
   int sigma_1[CODE_LENGTH_N];
   unsigned char seedPath[LOGN][SEEDBYTES];
   unsigned char seeds[N][SEEDBYTES];

   /* Parse response */
   unpack_response(resp, i, commitments[i], tilde_e_list[i],
                   sizeof(tilde_e_list[i]), sigma_1, sizeof(sigma_1),
                   seedPath);
   genSeeds(seeds, seedPath, i);

   if (i != 0) {
      /* Check if sigma1 is valid */
      if (1 != verify_G(sigma_1, non_sys_null_G)) {
         fprintf(stderr, "REJECT: sigma1 is not valid!\n");
         return 1;
      }
   }

   int v_list[N][CODE_LENGTH_N];
   int sigma_list[N][CODE_LENGTH_N];
   unsigned char salt_list[N][SEEDBYTES];

   //Save sigma1
   if (i != 0) {
      for (int j = 0; j < CODE_LENGTH_N; j++)
         sigma_list[0][j] = sigma_1[j];
   }

   for (int j = 0; j < N; j++) {
      if (j == i)
         continue;

      //For j = 0, we sample only v
      if (j == 0) {
         sample_from_seed_only_v(seeds[0], v_list[0], salt_list[0]);
      } else {
         sample_from_seed(seeds[j], non_sys_base_G, sigma_list[j],
                          v_list[j], salt_list[j]);
      }

#ifdef PRINTING
      printf("verif_salt{1,%d} = '", j + 1);
      for (int k = 0; k < SEEDBYTES; k++) {
         printf("%02x", salt_list[j][k]);
      }
      printf("%02x';", salt_list[j][SEEDBYTES - 1]);

      printf("verif_v_list{1,%d} = [", j + 1);
      for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
         printf("%d, ", v_list[j][k]);
      }
      printf("%d];", v_list[j][CODE_LENGTH_N - 1]);

      printf("verif_sigma_list{1,%d} = [", j + 1);
      for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
         printf("%d, ", sigma_list[j][k]);
      }
      printf("%d];", sigma_list[j][CODE_LENGTH_N - 1]);
#endif
   }

   //Recompute tilde_e_0 and tilde_e_1, only if i != 0
   int tmp_tilde_e[CODE_LENGTH_N];
   if (i != 0) {
      for (int i = 0; i < CODE_LENGTH_N; i++) {
         tmp_tilde_e[i] = (beta * square_and_multiply(e[i])) % Fq_size;
      }
      multiply_restricted_by_Fq(sigma_list[0], tmp_tilde_e, tmp_tilde_e);
      sum_Fq(tmp_tilde_e, v_list[0], tilde_e_list[0]);
   }
#ifdef PRINTING
   printf("verif_tilde_e0 = [");
   for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
      printf("%d, ", tmp_tilde_e[k]);
   }
   printf("%d];\n", tmp_tilde_e[CODE_LENGTH_N - 1]);

   printf("verif_tilde_e{1, 1} = [");
   for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
      printf("%d, ", tilde_e_list[0][k]);
   }
   printf("%d]\n;", tilde_e_list[0][CODE_LENGTH_N - 1]);
#endif

   /* ej */
   for (int j = 1; j < N; j++) {
      if (j == i)
         continue;
      multiply_restricted_by_Fq(sigma_list[j], tilde_e_list[j - 1],
                                tilde_e_list[j]);
      sum_Fq(tilde_e_list[j], v_list[j], tilde_e_list[j]);


#ifdef PRINTING
      printf("verif_tilde_e{1,%d} = [", j + 1);
      for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
         printf("%d, ", tilde_e_list[j][k]);
      }
      printf("%d];\n", tilde_e_list[j][CODE_LENGTH_N - 1]);
#endif
   }

#ifdef PRINTING
   printf("verif_tilde_e{1,%d} = [", i + 1);
   for (int k = 0; k < CODE_LENGTH_N - 1; k++) {
      printf("%d, ", tilde_e_list[i][k]);
   }
   printf("%d];\n", tilde_e_list[i][CODE_LENGTH_N - 1]);
#endif

   /* Regenerate commitments */
   unsigned char commit_input[2 * SEEDBYTES + CODE_LENGTH_N];
   for (int j = 0; j < N; j++) {
      if (j == i)
         continue;

      /* Now recompute commitment */
      memcpy(commit_input, salt_list[j], SEEDBYTES);
      memcpy(&commit_input[SEEDBYTES], seeds[j], SEEDBYTES);

      if (j != 0) {
         sha3_256(commitments[j],commit_input, 2 * SEEDBYTES );
      } else {

         /* The following piece of code is not efficient! */
         for (int k = 0; k < SEEDBYTES; k++)
            commit_input[k] = seeds[0][k];
         for (int k = 0; k < SEEDBYTES; k++)
            commit_input[k + SEEDBYTES] = salt_list[0][k];
         for (int k = 0; k < CODE_LENGTH_N; k++)
            commit_input[k + 2 * SEEDBYTES] = (unsigned char) sigma_list[0][k];

         /*
            memcpy(&commit_input[0], seeds[0], SEEDBYTES);
            memcpy(&commit_input[SEEDBYTES], salt_list[0], SEEDBYTES);
            memcpy(&commit_input[2*SEEDBYTES], sigma_1, CODE_LENGTH_N);
          */

#ifdef PRINTING
         fprintf(stderr, "verif_c1_input = '");
         for (int k = 0; k < 2 * SEEDBYTES + CODE_LENGTH_N - 1; k++)
            fprintf(stderr, "%02x", commit_input[k]);
         fprintf(stderr, "%02x';\n",
                 commit_input[2 * SEEDBYTES + CODE_LENGTH_N - 1]);
#endif

         sha3_256(commitments[0],commit_input, 2 * SEEDBYTES + CODE_LENGTH_N);

#ifdef PRINTING
         fprintf(stderr, "c{1,1} = '");
         for (int k = 0; k < HASHBYTES - 1; k++)
            fprintf(stderr, "%02x", commitments[0][k]);
         fprintf(stderr, "%02x';\n", commitments[0][HASHBYTES - 1]);
#endif
      }
   }

   for (int j = 0; j < N; j++) {
#ifdef PRINTING
      fprintf(stderr, "verif_c{1,%d} = '", j + 1);
      for (int k = 0; k < HASHBYTES - 1; k++)
         fprintf(stderr, "%02x", commitments[j][k]);
      fprintf(stderr, "%02x';\n", commitments[j][HASHBYTES - 1]);
#endif
   }
   /*      Compute syndrome */
   int ver_s[CODE_LENGTH_N - CODE_DIMENSION_K];

   /* sigma(e)*H^T */
   compute_syndrome_Fq(tilde_e_list[N - 1], non_sys_H, ver_s);

   /* sigma(e)*H^T - beta*s */
   for (int j = 0; j < CODE_LENGTH_N - CODE_DIMENSION_K; j++) {
      int ver_s_j = (ver_s[j] - beta * s[j]);
      while (ver_s_j < 0)
         ver_s_j += Fq_size;
      ver_s[j] = ver_s_j;
   }

#ifdef PRINTING
   fprintf(stderr, "verif_s = [");
   for (int k = 0; k < CODE_LENGTH_N - CODE_DIMENSION_K - 1; k++)
      fprintf(stderr, "%d, ", ver_s[k]);
   fprintf(stderr, "%d];\n", ver_s[CODE_LENGTH_N - CODE_DIMENSION_K - 1]);
#endif

   /* Hash syndrome and commitments */
   unsigned char c_input[2 * (CODE_LENGTH_N - CODE_DIMENSION_K) + N * HASHBYTES];

   /* Convert to char array */
   unsigned char ver_s_char[2 * (CODE_LENGTH_N - CODE_DIMENSION_K)];
   unsigned char ver_c[HASHBYTES];
   array2string(ver_s, CODE_LENGTH_N - CODE_DIMENSION_K, ver_s_char);

   memcpy(c_input, ver_s_char, 2 * (CODE_LENGTH_N - CODE_DIMENSION_K));
   for (int k = 0; k < N; k++) {
      memcpy(&c_input
             [2 * (CODE_LENGTH_N - CODE_DIMENSION_K) + k * HASHBYTES],
             &commitments[k], HASHBYTES);
   }
   sha3_256(ver_c,c_input, HASHBYTES);

#ifdef PRINTING
   fprintf(stderr, "verif_final_c = '");
   for (int k = 0; k < HASHBYTES - 1; k++)
      fprintf(stderr, "%02x", ver_c[k]);
   fprintf(stderr, "%02x';\n", ver_c[HASHBYTES - 1]);
#endif

   /* Hash tilde_e */
   unsigned char h_input[2 * N * CODE_LENGTH_N];

   /* Convert to char array */
   unsigned char tilde_e_char[2 * CODE_LENGTH_N];
   unsigned char ver_h[HASHBYTES];
   for (int k = 0; k < N; k++) {
      array2string(tilde_e_list[k], CODE_LENGTH_N, tilde_e_char);
//        memcpy(&h_input[k*CODE_LENGTH_N], tilde_e_char, 2*CODE_LENGTH_N);
      for (int u = 0; u < 2 * CODE_LENGTH_N; u++) {
         h_input[k * 2 * CODE_LENGTH_N + u] = tilde_e_char[u];
      }
   }

   sha3_256(ver_h,h_input, HASHBYTES);

#ifdef PRINTING
   fprintf(stderr, "verif_final_h = '");
   for (int k = 0; k < HASHBYTES - 1; k++)
      fprintf(stderr, "%02x", ver_h[k]);
   fprintf(stderr, "%02x';\n", ver_h[HASHBYTES - 1]);
#endif


   int ok = 0;
   for (int k = 0; k < HASHBYTES; k++) {
      ok |= (ver_c[k] ^ c[k]) | (ver_h[k] ^ h[k]);
   }

   /* returns 1 if ver_c and c or ver_h and h don't match */
   return ok;
}
