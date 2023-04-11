#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"
#include "utils.h"
#include "time.h"
#include "rng.h"
#include "string.h"

#include "keygen.h"
#include "commit.h"
#include "response.h"
#include "verify.h"
#include "keccak.h"

#define ROUNDS 10000

int main(int argc, char *argv[]) {
	
	FILE* fout; //file for printing results

	
	fprintf(stderr,"Considered parameters q = %d, g = %d, Ord(g) = z = %d, n = %d, k = %d, m = %d, N = %d \n", Fq_size, E_gen, Z, code_length, code_dimension, M, N);
	
	double cpu_time_used;
	
	//Variables for time measurements     
    clock_t t0, t1, t2, t3, t4, t5;
	int keygen_clocks = 0, commit_clocks = 0, secondcommit_clocks = 0, response_clocks = 0, verification_clocks = 0;
		
	//Useful to have reproducible results
    srand(0);

	for(int num=0; num<ROUNDS; num++){
		
		t0 = clock();
		
	    //Objects in key pair 
		int non_sys_H[code_length-code_dimension][code_dimension];
		int non_sys_base_G[M ][code_length - M], non_sys_null_G[code_length - M][M];
		int e[code_length], sigma[code_length], sigma_e[code_length]; 
		int s[code_length - code_dimension]; 
		
		//Generate key pair	
		key_gen(non_sys_H, non_sys_base_G, non_sys_null_G, e, sigma, sigma_e, s);
		
		t1 = clock();
		
		//Uncomment if you want to print key pair
		#ifdef PRINTING
			print_key_pair(non_sys_H, non_sys_base_G, non_sys_null_G, e, sigma, sigma_e, s);
		#endif
			
	    
		//Start commitments
		int v_list[N][code_length]; //matrix with vector v_1, v_2, ..., v_N
		int sigma_list[N][code_length]; //matrix with transformations sigma_1, sigma_2, ..., sigma_N
		char commitments[N][HASHBYTES];
		char seedTree[2*N-1][SEEDBYTES], c[HASHBYTES], h[HASHBYTES];
		int accumulate_sigma[code_length]={0}; //useful for challenge phase

        /* First commitment phase */
		commit(sigma, non_sys_H, non_sys_base_G, v_list, sigma_list, commitments, accumulate_sigma, c, seedTree);
		
		t2 = clock();
		
		//Uncomment if you want to print v and sigma
		#ifdef PRINTING
			print_v_and_sigma(v_list, sigma_list);
		#endif		
		/* Sample beta */
		int beta = 1+(rand()%(Fq_size-1));
		
		/* Second commitment phase */
		int tilde_e_list[N][code_length];
			
        /* Added h as hash output here and in second_commit. */
        /* Could actually reuse c as buffer, but for clarity this way is more consistent with the pseudocode */
		second_commit(beta, e, sigma, non_sys_H, non_sys_base_G, v_list, sigma_list, accumulate_sigma, tilde_e_list, h);
		
		t3 = clock();
		
        //Uncomment if you want to print e_tilde
        /* print_e_tilde(tilde_e_list); */
       	#ifdef PRINTING
 	       print_e_tilde(tilde_e_list);	
		#endif
		
        /* Sample challenge */
        int i = rand()%N;
        /* printf("i= %d\n", i); */

        /* Response */
        Resp resp;
        char seed_path[LOGN][SEEDBYTES];
        gen_seed_path(seed_path, seedTree, i);
        pack_response(&resp, i, commitments[i], tilde_e_list[i], sizeof(tilde_e_list[i]), sigma_list[0], sizeof(sigma_list[0]), seed_path); 
		
		t4 = clock();
		
        /* Verification */
        if ( verify( &resp, i, non_sys_null_G, non_sys_base_G, non_sys_H, beta, e, s, c, h) ) {
            fprintf(stderr,"REJECT %d: verification failed - c and/or h do not match!\n", num);
        } else {
            //fprintf(stderr,"ACCEPT %d: verification successful!\n", num);
        }
        
        t5 = clock();

		keygen_clocks += (t1 - t0);
		commit_clocks += (t2 - t1);
		secondcommit_clocks += (t3 - t2);
		response_clocks += (t4 - t3);
		verification_clocks += (t5 - t4);
		
	}
	

	fout = fopen("output.txt", "w");
	fprintf(fout, "Number of tests: %d\n", ROUNDS);
	fprintf(fout, "Key generation: average number of clock cycles (kCycles)) = %f, ms = %f\n", (double) keygen_clocks / ROUNDS, (double) 1000* keygen_clocks / CLOCKS_PER_SEC / ROUNDS);
	fprintf(fout, "Commitment: average number of clock cycles (kCycles)) = %f, ms = %f\n", (double) commit_clocks / ROUNDS, (double)1000 * commit_clocks / CLOCKS_PER_SEC / ROUNDS);
	fprintf(fout, "Second commitment: average number of clock cycles (kCycles)) = %f, ms = %f\n", (double) secondcommit_clocks / ROUNDS, (double)1000 * secondcommit_clocks / CLOCKS_PER_SEC / ROUNDS);
	fprintf(fout, "Response: average number of clock cycles (kCycles)) = %f, ms = %f\n", (double) response_clocks / ROUNDS, (double)1000 * response_clocks / CLOCKS_PER_SEC / ROUNDS);
	fprintf(fout, "Verification: average number of clock cycles (kCycles)) = %f, ms = %f\n", (double) verification_clocks / ROUNDS, (double)1000 * verification_clocks / CLOCKS_PER_SEC / ROUNDS);
	fclose(fout);

	system("PAUSE");
		      	
    return 0;
}

