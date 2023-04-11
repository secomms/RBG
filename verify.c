#include "verify.h"

int verify( const Resp *resp, const int i, const int non_sys_null_G[code_length - M][M], const int non_sys_base_G[M][code_length - M], 
            const int non_sys_H[code_length-code_dimension][code_dimension], const int beta, const int e[code_length], const int s[code_length - code_dimension],
            const char c[HASHBYTES], const char h[HASHBYTES])
{
    char commitments[N][HASHBYTES];
    int tilde_e_list[N][code_length], tilde_e_0[code_length];
    int sigma_1[code_length];
    char seedPath[LOGN][SEEDBYTES];
    char seeds[N][SEEDBYTES];

    /* Parse response */
    unpack_response(resp, i, commitments[i], tilde_e_list[i], sizeof(tilde_e_list[i]), sigma_1, sizeof(sigma_1), seedPath);
    genSeeds(seeds, seedPath, i);

    if (i != 0) {
        /* Check if sigma1 is valid */
        if (1 != verify_G(sigma_1, non_sys_null_G)) { 
            fprintf(stderr,"REJECT: sigma1 is not valid!\n");
            return 1;
        }
    }
		
    int v_list[N][code_length];
    int sigma_list[N][code_length];
    char salt_list[N][SEEDBYTES];

	//Save sigma1
	if(i!=0){
        for(int j = 0;j<code_length;j++)sigma_list[0][j] = sigma_1[j];
	}
	
    for (int j=0; j<N; j++) {
        if (j == i) continue;
        
        //For j = 0, we sample only v
        if(j==0){
        	sample_from_seed_only_v(seeds[0], v_list[0], salt_list[0]);
		}else{
			sample_from_seed(seeds[j], non_sys_base_G, sigma_list[j], v_list[j], salt_list[j]);
		}
		
    	#ifdef PRINTING			
	        printf("verif_salt{1,%d} = '",j+1);
	        for (int k=0; k<SEEDBYTES; k++){ 
				printf("%02x", salt_list[j][k]); 
			}
			printf("%02x';",salt_list[j][SEEDBYTES-1]); 			
			
			printf("verif_v_list{1,%d} = [",j+1);
	        for (int k=0; k<code_length-1; k++) { 
				printf("%d, ", v_list[j][k]); 
			}
			printf("%d];", v_list[j][code_length-1]); 			
			
			printf("verif_sigma_list{1,%d} = [",j+1);
	        for (int k=0; k<code_length-1; k++) { 
				printf("%d, ", sigma_list[j][k]); 
			}
			printf("%d];", sigma_list[j][code_length-1]); 			
		#endif
    }
	
	//Recompute tilde_e_0 and tilde_e_1, only if i != 0
	int tmp_tilde_e[code_length];
	if(i!=0){
		for(int i = 0;i<code_length;i++){
			tmp_tilde_e[i] = (beta*square_and_multiply(e[i]))%Fq_size;
		}
		multiply_restricted_by_Fq(sigma_list[0], tmp_tilde_e, tmp_tilde_e);
		sum_Fq(tmp_tilde_e, v_list[0], tilde_e_list[0]);
	}
	
	#ifdef PRINTING			
		printf("verif_tilde_e0 = [");
	    for (int k=0; k<code_length-1; k++){
			printf("%d, ", tmp_tilde_e[k]); 
		}
		printf("%d];\n", tmp_tilde_e[code_length-1]); 			
			
		printf("verif_tilde_e{1, 1} = [");
	    for (int k=0; k<code_length-1; k++){
			printf("%d, ", tilde_e_list[0][k]); 
		}
		printf("%d]\n;", tilde_e_list[0][code_length-1]);
	#endif
		
    /* ej */
    for (int j=1; j<N; j++){
        if (j == i) continue;
        multiply_restricted_by_Fq(sigma_list[j], tilde_e_list[j-1], tilde_e_list[j]);
        sum_Fq(tilde_e_list[j], v_list[j], tilde_e_list[j]);
        
        
        #ifdef PRINTING			
	        printf("verif_tilde_e{1,%d} = [",j+1);
		    for (int k=0; k<code_length-1; k++){
				printf("%d, ", tilde_e_list[j][k]); 
			}
			printf("%d];\n", tilde_e_list[j][code_length-1]);
		#endif
    }
    
    #ifdef PRINTING			
	    printf("verif_tilde_e{1,%d} = [",i+1);
	    for (int k=0; k<code_length-1; k++){
			printf("%d, ", tilde_e_list[i][k]); 
		}
		printf("%d];\n", tilde_e_list[i][code_length-1]);
	#endif
			
    /* Regenerate commitments */
    char commit_input[2*SEEDBYTES+code_length];
    for (int j=0; j<N; j++) {
        if (j == i) continue;

        /* Now recompute commitment */
        memcpy(commit_input, salt_list[j], SEEDBYTES);
        memcpy(&commit_input[SEEDBYTES], seeds[j], SEEDBYTES);

        if (j != 0){
            FIPS202_SHA3_256(commit_input, 2*SEEDBYTES, commitments[j]);
        } else {

			/* The following piece of code is not efficient!*/
			for(int k = 0;k<SEEDBYTES;k++)commit_input[k] = seeds[0][k];
			for(int k = 0;k<SEEDBYTES;k++)commit_input[k+SEEDBYTES] = salt_list[0][k];
			for(int k = 0;k<code_length;k++)commit_input[k+2*SEEDBYTES] = (char)sigma_list[0][k];
		
			/*
        	memcpy(&commit_input[0], seeds[0], SEEDBYTES);
        	memcpy(&commit_input[SEEDBYTES], salt_list[0], SEEDBYTES);
            memcpy(&commit_input[2*SEEDBYTES], sigma_1, code_length);
            */
            
        	#ifdef PRINTING			
	            fprintf(stderr,"verif_c1_input = '");
				for(int k = 0;k<2*SEEDBYTES+code_length-1;k++)fprintf(stderr,"%02x",commit_input[k]);
				fprintf(stderr,"%02x';\n",commit_input[2*SEEDBYTES+code_length-1]);
			#endif
	
            FIPS202_SHA3_256(commit_input, 2*SEEDBYTES+code_length, commitments[0]);
            
        	#ifdef PRINTING			
	            fprintf(stderr,"c{1,1} = '");
				for(int k = 0;k<HASHBYTES-1;k++)fprintf(stderr,"%02x",commitments[0][k]);
				fprintf(stderr,"%02x';\n",commitments[0][HASHBYTES-1]);
			#endif
        }
    }

	for(int j = 0;j<N;j++){
		#ifdef PRINTING			
			fprintf(stderr,"verif_c{1,%d} = '",j+1);
			for(int k = 0;k<HASHBYTES-1;k++)fprintf(stderr,"%02x",commitments[j][k]);
			fprintf(stderr,"%02x';\n",commitments[j][HASHBYTES-1]);
		#endif
	}
	/* 	Compute syndrome */	
    int ver_s[code_length-code_dimension];

    /* sigma(e)*H^T */
    compute_syndrome_Fq(tilde_e_list[N-1], non_sys_H, ver_s);

    /* sigma(e)*H^T - beta*s*/
    for(int j = 0;j<code_length-code_dimension;j++){
        int ver_s_j = (ver_s[j] - beta*s[j]);
        while(ver_s_j<0) ver_s_j+=Fq_size;
        ver_s[j] = ver_s_j;
    }
	
	#ifdef PRINTING			
		fprintf(stderr,"verif_s = [");
		for(int k = 0;k<code_length - code_dimension -1;k++)fprintf(stderr,"%d, ",ver_s[k]);
		fprintf(stderr,"%d];\n",ver_s[code_length - code_dimension-1]);
	#endif
	
    /* Hash syndrome and commitments */
    char c_input[2*(code_length-code_dimension)+N*HASHBYTES];

    /* Convert to char array */
    char ver_s_char[2*(code_length-code_dimension)];
    char ver_c[HASHBYTES];
    array2string(ver_s, code_length-code_dimension, ver_s_char);

    memcpy(c_input, ver_s_char, 2*(code_length-code_dimension));
    for (int k=0; k<N; k++){
        memcpy(&c_input[2*(code_length-code_dimension)+k*HASHBYTES], &commitments[k], HASHBYTES);
    }
	FIPS202_SHA3_256(c_input, HASHBYTES, ver_c);
	
	#ifdef PRINTING			
		fprintf(stderr,"verif_final_c = '");
		for(int k = 0;k<HASHBYTES-1;k++)fprintf(stderr,"%02x",ver_c[k]);
		fprintf(stderr,"%02x';\n",ver_c[HASHBYTES-1]);
	#endif
		
    /* Hash tilde_e */
    char h_input[2*N*code_length];

    /* Convert to char array */
    char tilde_e_char[2*code_length];
    char ver_h[HASHBYTES];
    for (int k=0; k<N; k++){
        array2string(tilde_e_list[k], code_length, tilde_e_char);
//        memcpy(&h_input[k*code_length], tilde_e_char, 2*code_length);
		for(int u = 0;u<2*code_length;u++){
			h_input[k*2*code_length + u] = tilde_e_char[u];		
		}
    }
    
	FIPS202_SHA3_256(h_input, HASHBYTES, ver_h);	
	
	#ifdef PRINTING			
		fprintf(stderr,"verif_final_h = '");
		for(int k = 0;k<HASHBYTES-1;k++)fprintf(stderr,"%02x",ver_h[k]);
		fprintf(stderr,"%02x';\n",ver_h[HASHBYTES-1]);
	#endif
	

    int ok = 0;
    for (int k=0; k<N; k++) {
        ok |= (ver_c[k] ^ c[k]) | (ver_h[k] ^ h[k]);        
    }

    /* returns 1 if ver_c and c or ver_h and h don't match */
    return ok;
}
