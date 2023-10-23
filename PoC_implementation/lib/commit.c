#include "commit.h"
#include "sha3.h"

void commit(const int sigma[CODE_LENGTH_N],
        const int non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K],
        const int non_sys_base_G[M][CODE_LENGTH_N - M],
        int v_list[N][CODE_LENGTH_N], int sigma_list[N][CODE_LENGTH_N],
        unsigned char commitments[N][2 * SEEDBYTES],
        int accumulate_sigma[CODE_LENGTH_N], 
        unsigned char c[HASHBYTES],
        unsigned char seedTree[2 * N - 1][SEEDBYTES]){

   /* sample mseed */
   sample_seed(seedTree[0]);
   gen_seed_tree(seedTree);


   /* Leafs of the seedTree start at position seedTree[N-1] */
   /* printf("Salts commit\n"); */
   /* printf("V commit\n"); */
   for (int i = 1; i < N; i++) {

      //Sample seed_i and salt_i
      unsigned char salt_i[SEEDBYTES];
      sample_from_seed(seedTree[N - 1 + i], non_sys_base_G,
                       sigma_list[i], v_list[i], salt_i);


#ifdef PRINTING
      printf("prov_salt{1,%d} = '", i + 1);
      for (int k = 0; k < SEEDBYTES; k++) {
         printf("%02x", salt_i[k]);
      }
      printf("%02x';", salt_i[SEEDBYTES - 1]);
#endif

      //Hash for commitment
      unsigned char commit_input[HASHBYTES];    //contains input for hash
      for (int u = 0; u < SEEDBYTES; u++) {
         commit_input[u] = salt_i[u];
         commit_input[SEEDBYTES + u] = seedTree[N - 1 + i][u];
      }
      sha3_256(commitments[i],commit_input,HASHBYTES);

#ifdef PRINTING
      fprintf(stderr, "c{1,%d} = '", i + 1);
      for (int j = 0; j < HASHBYTES - 1; j++)
         fprintf(stderr, "%02x", commitments[i][j]);
      fprintf(stderr, "%02x';\n", commitments[i][HASHBYTES - 1]);
#endif

      //Accumulate restricted isometries (useful later)
      multiply_restricted(accumulate_sigma, sigma_list[i],
                          accumulate_sigma);

   }


   //Sample v_1 and sigma_1
   unsigned char salt_1[SEEDBYTES];

   int v[CODE_LENGTH_N];

   sample_from_seed_only_v(seedTree[N - 1], v_list[0],
                           salt_1);    //sample only v and salt
   divide_restricted(sigma, accumulate_sigma, sigma_list[0]);    //compute sigma_1

   //Compute c1
   unsigned char c1_input[2 * SEEDBYTES + CODE_LENGTH_N];
   for (int u = 0; u < SEEDBYTES; u++) {
      c1_input[u] = seedTree[N - 1][u];
      c1_input[u + SEEDBYTES] = salt_1[u];
   }

   for (int u = 0; u < CODE_LENGTH_N; u++) {
      c1_input[2 * SEEDBYTES + u] = (unsigned char) (sigma_list[0][u]);
   }


   sha3_256(commitments[0], c1_input, HASHBYTES);
#ifdef PRINTING
   fprintf(stderr, "c{1,1} = '");
   for (int k = 0; k < 2 * SEEDBYTES + CODE_LENGTH_N - 1; k++)
      fprintf(stderr, "%02x", commitments[0][k]);
   fprintf(stderr, "%02x';\n",
           commitments[0][2 * SEEDBYTES + CODE_LENGTH_N - 1]);
#endif

   //Start computing v
   multiply_restricted_by_Fq(accumulate_sigma, v_list[0], v);

   //Compute v
   for (int i = 1; i < N - 1; i++) {
      divide_restricted(accumulate_sigma, sigma_list[i],
                        accumulate_sigma);    //divide by sigma_{i+1}
      int this_v[CODE_LENGTH_N];
      multiply_restricted_by_Fq(accumulate_sigma, v_list[i], this_v);
      sum_Fq(v, this_v, v);
   }
   sum_Fq(v, v_list[N - 1], v);    //sum v_N

   int sin_v[CODE_LENGTH_N - CODE_DIMENSION_K];
   compute_syndrome_Fq(v, non_sys_H, sin_v);
#ifdef PRINTING
   fprintf(stderr, "sin_v = [");
   for (int j = 0; j < CODE_LENGTH_N - CODE_DIMENSION_K - 1; j++)
      fprintf(stderr, "%d, ", sin_v[j]);
   fprintf(stderr, "%d];\n", sin_v[CODE_LENGTH_N - CODE_DIMENSION_K - 1]);
#endif

   //Do final hash
   unsigned char c_input[2 * (CODE_LENGTH_N - CODE_DIMENSION_K) + N * HASHBYTES];

   unsigned char sin_v_char[2 * (CODE_LENGTH_N - CODE_DIMENSION_K)];
   array2string(sin_v, CODE_LENGTH_N - CODE_DIMENSION_K, sin_v_char);
   for (int u = 0; u < 2 * (CODE_LENGTH_N - CODE_DIMENSION_K); u++) {
      c_input[u] = sin_v_char[u];
   }

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < HASHBYTES; j++) {
         c_input[2 * (CODE_LENGTH_N - CODE_DIMENSION_K) + i * HASHBYTES +
                 j] = commitments[i][j];
      }
   }

   sha3_256(c,c_input, HASHBYTES);

#ifdef PRINTING
   fprintf(stderr, "final_c = '");
   for (int k = 0; k < HASHBYTES - 1; k++)
      fprintf(stderr, "%02x", c[k]);
   fprintf(stderr, "%02x';\n", c[HASHBYTES - 1]);
#endif

   return;
}






void print_v_and_sigma(int v_list[N][CODE_LENGTH_N],
                       int sigma_list[N][CODE_LENGTH_N])
{


   fprintf(stderr, "v_list{1,%d} = [", 1);
   for (int j = 0; j < CODE_LENGTH_N - 1; j++)
      fprintf(stderr, "%d, ", v_list[0][j]);
   fprintf(stderr, "%d];\n", v_list[0][CODE_LENGTH_N - 1]);

   fprintf(stderr, "sigma_list{1,%d} = [", 1);
   for (int j = 0; j < CODE_LENGTH_N - 1; j++)
      fprintf(stderr, "%d, ", sigma_list[0][j]);
   fprintf(stderr, "%d];\n", sigma_list[0][CODE_LENGTH_N - 1]);

   for (int i = 1; i < N; i++) {

      fprintf(stderr, "v_list{1,%d} = [", i + 1);
      for (int j = 0; j < CODE_LENGTH_N - 1; j++)
         fprintf(stderr, "%d, ", v_list[i][j]);
      fprintf(stderr, "%d];\n", v_list[i][CODE_LENGTH_N - 1]);

      fprintf(stderr, "sigma_list{1,%d} = [", i + 1);
      for (int j = 0; j < CODE_LENGTH_N - 1; j++)
         fprintf(stderr, "%d, ", sigma_list[i][j]);
      fprintf(stderr, "%d];\n", sigma_list[i][CODE_LENGTH_N - 1]);

   }

   return;
}

void second_commit(const int beta, const int e[CODE_LENGTH_N],
                   const int sigma[CODE_LENGTH_N],
                   const int non_sys_H[CODE_LENGTH_N -
                                       CODE_DIMENSION_K][CODE_DIMENSION_K],
                   const int non_sys_base_G[M][CODE_LENGTH_N - M],
                   const int v_list[N][CODE_LENGTH_N],
                   const int sigma_list[N][CODE_LENGTH_N],
                   int accumulate_sigma[CODE_LENGTH_N],
                   int tilde_e_list[N][CODE_LENGTH_N], 
                   unsigned char h[HASHBYTES])
{

   //Start doing tilde_e_i
   int tmp_tilde_e[CODE_LENGTH_N];
   for (int i = 0; i < CODE_LENGTH_N; i++) {
      tmp_tilde_e[i] = (beta * square_and_multiply(e[i])) % Fq_size;
   }



   multiply_restricted_by_Fq(sigma_list[0], tmp_tilde_e, tmp_tilde_e);
   sum_Fq(tmp_tilde_e, v_list[0], tilde_e_list[0]);

   for (int i = 1; i < N; i++) {

      multiply_restricted_by_Fq(sigma_list[i], tilde_e_list[i - 1],
                                tmp_tilde_e);    //Apply sigma_i
      sum_Fq(tmp_tilde_e, v_list[i], tilde_e_list[i]);    //sum v_i

   }

   unsigned char h_input[N * 2 * CODE_LENGTH_N];

   for (int i = 0; i < N; i++) {
      unsigned char tilde_e_string[2 * CODE_LENGTH_N];
      array2string(tilde_e_list[i], CODE_LENGTH_N, tilde_e_string);
      for (int j = 0; j < CODE_LENGTH_N; j++) {
         h_input[i * 2 * CODE_LENGTH_N + j] = tilde_e_string[j];
      }
   }
   sha3_256(h, h_input, HASHBYTES);

#ifdef PRINTING
   fprintf(stderr, "final_h = '");
   for (int k = 0; k < HASHBYTES - 1; k++)
      fprintf(stderr, "%02x", h[k]);
   fprintf(stderr, "%02x';\n", h[HASHBYTES - 1]);
#endif

   return;
}

void print_e_tilde(int tilde_e_list[N][CODE_LENGTH_N])
{
   for (int i = 0; i < N; i++) {
      fprintf(stderr, "tilde_e{1,%d} = [", i + 1);
      for (int j = 0; j < CODE_LENGTH_N - 1; j++) {
         fprintf(stderr, "%d, ", tilde_e_list[i][j]);
      }
      fprintf(stderr, "%d];\n", tilde_e_list[i][CODE_LENGTH_N - 1]);
   }
   return;
}
