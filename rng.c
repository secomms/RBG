#include "rng.h"

//generate a random seed and salt
void sample_seed_and_salt (char stream[2*SEEDBYTES])
{

  for (int i = 0; i < 2*SEEDBYTES; i++)
  {
  	int x = rand();
    stream[i] = x;
  }

  return ;
}


/*********************************************/

//generate a random seed
void sample_seed (char stream[SEEDBYTES])
{

  for (int i = 0; i < SEEDBYTES; i++)
  {
  	int x = rand();
    stream[i] = x;
  }

  return ;
}


/*************************************************************/


/**************************************************/

void sample_from_seed(const char seed[SEEDBYTES], const int non_sys_base_G[M][code_length - M], int sigma[code_length], int v[code_length], char salt[SEEDBYTES]){
	
	
	char digest[NUMBYTES+SEEDBYTES];
	FIPS202_SHAKE128(seed, SEEDBYTES, digest, NUMBYTES+SEEDBYTES);    


	//parse digest to create sigma
	for(int i = 0;i<code_length-M;i++)sigma[M+i] = 0;
	
    for(int i = 0; i < M; i++){
    	int val = 0;
    	char ch = digest[i];
    	for(int j=0;j<8;j++){
    		int x = ch&(1<<j);
    		if(x){
    			val+=(1<<j);
			}else{
				val+=0;
			}
		}
		val = val%Z;
		
		//Update values of sigma
		sigma[i] = val;		
    	for (int j = 0;j<code_length-M;j++){
    		sigma[M+j]=(sigma[M+j]+val*non_sys_base_G[i][j])%Z;
		}
    }
	
	
	//parse digest to create v
    for(int i = 0; i < code_length; i++){
    	int val = 0;
    	char ch = digest[M+2*i];
    	for(int j=0;j<8;j++){
    		int x = ch&(1<<j);
    		if(x){
    			val+=(1<<j);
			}else{
				val+=0;
			}
		}

		ch = digest[M+2*i+1];
    	for(int j=0;j<8;j++){
    //	for(int j=0;j<2;j++){
    		int x = ch&(1<<(j));
    		if(x){
    			val+=(1<<(8+j));
			}else{
				val+=0;
			}
		}
		v[i] = val%Fq_size;
    }
	
	for(int i = 0;i<SEEDBYTES;i++)salt[i] = digest[M+2*code_length+i];
	
	return;	
}
	

/***************************************************************/
void sample_from_seed_only_v(const char seed[SEEDBYTES], int v[code_length], char salt[SEEDBYTES]){
	
	
    char digest[2*code_length+SEEDBYTES];
	FIPS202_SHAKE128(seed, SEEDBYTES, digest, 2*code_length+SEEDBYTES);    
    /* char digest[NUMBYTES+SEEDBYTES]; */
	/* FIPS202_SHAKE128(seed, SEEDBYTES, digest, NUMBYTES+SEEDBYTES); */    
	
	//parse digest to create v
    for(int i = 0; i < code_length; i++){
    	int val = 0;
    	char ch = digest[2*i];
    	for(int j=0;j<8;j++){
    		int x = ch&(1<<j);
    		if(x){
    			val+=(1<<j);
			}else{
				val+=0;
			}
		}

		ch = digest[2*i+1];
    	for(int j=0;j<8;j++){
  //  	for(int j=0;j<2;j++){
    		int x = ch&(1<<(j));
    		if(x){
    			val+=(1<<(8+j));
			}else{
				val+=0;
			}
		}
		v[i] = val%Fq_size;
    }
	
	for(int i = 0;i<SEEDBYTES;i++)salt[i] = digest[2*code_length+i];
	/* for(int i = 0;i<SEEDBYTES;i++)salt[i] = digest[M+2*code_length+i]; */
		
	return;	
}

/* Generate the full seed tree (as in BG) including internal nodes (will be required for response, idc about memory consumption for now) */
/* Requires 2N-1 * SEEDBYTES bytes */
/***************************************************************/
void gen_seed_tree(char seedTree[2*N-1][SEEDBYTES]) 
{
    /* seedTree[0] already initialized with mseed */
    /* TODO: add prefix h to input, do this when using incremental API. Shouldn't be relevant for performance at the moment as it's only 1 byte and */ 
    /* no additional permutation is required */
    char tmp[2*SEEDBYTES];
    for (int u=0; u<LOGN; u++) {
        for (int l=0; l<(1<<u); l++) {
		    FIPS202_SHA3_256(seedTree[((1<<u)+l-1)], SEEDBYTES, tmp);	
            memcpy(seedTree[((1<<(u+1))+2*l-1)], tmp, SEEDBYTES);
            memcpy(seedTree[((1<<(u+1))+2*l)], &tmp[SEEDBYTES], SEEDBYTES);
        }
    }
    /* printf("\nSampled seedTree\n"); */
    /* for (int x=N-1; x<2*N-1; x++) { */
    /*     for (int y=0; y<SEEDBYTES; y++) { */
    /*         printf("%02x", (char)seedTree[x][y]); */
    /*     } */
    /*     printf("\n"); */
    /* } */
}

void gen_seed_path(char seedPath[LOGN][SEEDBYTES], const char seedTree[2*N-1][SEEDBYTES], const int i)
{
    int n;
    int node;
    int offset;

    /* Formula to compute path node seedPath[k] = n - 1 + i + offset */
    /* with the following initial values: */
    n       = N;
    node    = i;
    offset = (node&1) ? -1 : 1;

    /* node starts from leaf_i and traverses the tree upwards */
    for (int k=0; k<LOGN; k++) {
        memcpy(seedPath[k], seedTree[n-1+node+offset], SEEDBYTES);
        n >>= 1;
        node >>= 1;
        offset = (node&1) ? -1 : 1;
    }
}

void genSeeds(char seeds[N][SEEDBYTES], const char seedPath[LOGN][SEEDBYTES], const int i)
{
    
    int n;
    int node;
    int offset;
    char halfTree[N-1][SEEDBYTES] = {0};
    char tmp[2*SEEDBYTES];
    int full_tree_pos, flattened_pos;

    /* Equation to compute path node seedPath[k] = n - 1 + i + offset */
    /* with the following initial values: */
    n       = N;
    node    = i;
    offset = (node&1) ? -1 : 1;

    /* Build the partial trees from the seeds and copy their leaves into seeds[N][SEEDBYTES] */
    for (int k=0; k<LOGN; k++) {
        memcpy(halfTree[0], seedPath[k], SEEDBYTES);

        /* Build the small tree */
        for (int u=0; u<k; u++) {
            for (int l=0; l<(1<<u); l++) {
                FIPS202_SHA3_256(halfTree[((1<<u)+l-1)], SEEDBYTES, tmp);	
                memcpy(halfTree[((1<<(u+1))+2*l-1)], tmp, SEEDBYTES);
                memcpy(halfTree[((1<<(u+1))+2*l)], &tmp[SEEDBYTES], SEEDBYTES);
            }
        }

        /* position in original full seed tree and then in flattened one for seeds */
        full_tree_pos = n - 1 + node + offset;
        flattened_pos = (1<<k)*(full_tree_pos+1)-N;

        /* Now copy into seeds */
        for (int m=0; m<(1<<k); m++) {
            memcpy(seeds[flattened_pos+m], halfTree[(1<<k)-1+m], SEEDBYTES);
        }

        n >>= 1;
        node >>= 1;
        offset = (node&1) ? -1 : 1;
    }

    /* printf("\nGenerated seeds\n"); */
    /* for (int x=0; x<N; x++) { */
    /*     for (int y=0; y<SEEDBYTES; y++) { */
    /*         printf("%02x", (char)seeds[x][y]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

}













