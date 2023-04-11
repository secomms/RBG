#include "keygen.h"

void key_gen(int non_sys_H[code_length-code_dimension][code_dimension], int non_sys_base_G[M ][code_length - M], int non_sys_null_G[code_length - M][M], int e[code_length], int sigma[code_length], int sigma_e[code_length], int s[code_length - code_dimension]){
	
	gen_non_sys_H(non_sys_H); //H in systematic form
	gen_non_sys_base_G(non_sys_base_G, non_sys_null_G); //representation for G (in systematic form)
	gen_restricted(e, non_sys_base_G); //restricted e
	gen_restricted(sigma, non_sys_base_G); //restricted sigma
	multiply_restricted(e, sigma, sigma_e); //sigma(e)
	compute_syndrome(sigma_e, non_sys_H, s); //sigma(e)*H^T
	
	return ;
}


void print_key_pair(int non_sys_H[code_length-code_dimension][code_dimension], int non_sys_base_G[M ][code_length - M], int non_sys_null_G[code_length - M][M], int e[code_length], int sigma[code_length], int sigma_e[code_length], int s[code_length - code_dimension]){

	fprintf(stderr,"\n non_sys_H = [");
	for(int i=0;i<code_length - code_dimension-1;i++){
		for(int j = 0;j<code_dimension-1;j++ )fprintf(stderr,"%d, ",non_sys_H[i][j]);
		fprintf(stderr,"%d;\n",non_sys_H[i][code_dimension-1]);
	}
	for(int j = 0;j<code_dimension - 1;j++ )fprintf(stderr,"%d, ",non_sys_H[code_length - code_dimension-1][j]);
	fprintf(stderr,"%d]\n",non_sys_H[code_length-code_dimension-1][code_dimension-1]);
	
	
	fprintf(stderr,"\n non_sys_G = [");
	for(int i=0;i<M-1;i++){
		for(int j = 0;j<code_length-M-1;j++ )fprintf(stderr,"%d, ",non_sys_base_G[i][j]);
		fprintf(stderr,"%d;\n",non_sys_base_G[i][code_length-M-1]);
	}
	for(int j = 0;j<code_length-M-1;j++ )fprintf(stderr,"%d, ",non_sys_base_G[M-1][j]);
	fprintf(stderr,"%d]\n",non_sys_base_G[M-1][code_length-M-1]);



	fprintf(stderr,"\n non_sys_null_G = [");
	for(int i=0;i<code_length - M-1;i++){
		for(int j = 0;j<M-1;j++ )fprintf(stderr,"%d, ",non_sys_null_G[i][j]);
		fprintf(stderr,"%d;\n",non_sys_null_G[i][M-1]);
	}
	for(int j = 0;j<M-1;j++ )fprintf(stderr,"%d, ",non_sys_null_G[code_length - M-1][j]);
	fprintf(stderr,"%d]\n",non_sys_null_G[code_length-M-1][M-1]);
	//end of print facilities
	
	
	
	///PRINTING
	fprintf(stderr,"e = [");
	for(int j = 0;j<code_length-1;j++ )fprintf(stderr,"%d, ", e[j]);
	fprintf(stderr,"%d];\n",e[code_length-1]);
	
	fprintf(stderr,"sigma = [");
	for(int j = 0;j<code_length-1;j++ )fprintf(stderr,"%d, ", sigma[j]);
	fprintf(stderr,"%d];\n",sigma[code_length-1]);
	
	fprintf(stderr,"sigma_e = [");
	for(int j = 0;j<code_length-1;j++ )fprintf(stderr,"%d, ", sigma_e[j]);
	fprintf(stderr,"%d];\n",sigma_e[code_length-1]);
	
	fprintf(stderr,"s = [");
	for(int j = 0;j<code_length-code_dimension - 1;j++ )fprintf(stderr,"%d, ", s[j]);
	fprintf(stderr,"%d];\n",s[code_length-code_dimension-1]);
	/////////
	
}
