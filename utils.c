#include "utils.h"

//Generate non systematic part of H
void gen_non_sys_H(int  non_sys_H[code_length-code_dimension][code_dimension]){
	
	for(int i = 0; i<(code_length-code_dimension);i++){
		for(int j = 0; j<code_dimension;j++){
			int val = rand() % Fq_size;
			non_sys_H[i][j] = val;
		}	
	}
	
	return;
}


/********************************************************/

//Generate non systematic part for base_of_G
void gen_non_sys_base_G(int  non_sys_base_G[M ][code_length - M], int  non_sys_null_G[code_length - M][M]){
	
	for(int i = 0; i<M;i++){
		for(int j = 0; j<(code_length - M);j++){
			int val = rand() % Z;
			non_sys_base_G[i][j] = val;
			non_sys_null_G[j][i] = Z -val;
		}	
	}
	
	return;
}

/*********************************************************/

//Verify a restricted object is in G
int verify_G(int x[code_length], const int non_sys_null_G[code_length - M][M]){
	
	int s[code_length - M];
	
	for(int i = 0;i<code_length-M; i++){
		int x_i = x[M+i]; 
		s[i] = x_i;
	}
	
	
    for(int i = 0;i<M;i++){
    	
		int x_i = x[i]; //from exp to standard representation    	
    	for (int j = 0;j<(code_length - M);j++){
    		s[j]+=(x_i*non_sys_null_G[j][i])%Z;
		}
	}
	
	for(int i = 0;i<code_length - M;i++)s[i] = s[i]%Z;
    
/*	fprintf(stderr,"\n s = [");
    for(int i = 0;i<code_length-M-1;i++)fprintf(stderr,"%d, ",s[i]);
    fprintf(stderr,"%d]\n",s[code_length-M-1]);*/
    
	int ok = 1;
	for(int i = 0;i<code_length-M;i++){
		if(s[i])ok = 0;
	}
	
    return ok;
}

/*********************************************************************/

void gen_Fq(int x[code_length], const int seed){
    
    srand(seed);
    for(int i = 0;i<code_length;i++){
    	int coeff = rand()%Fq_size;
    	x[i] = coeff;
	}
	return ;
}

/********************************************************/

//Generate restricted object
//The returned vector will be a series of exponents to assign to g
void gen_restricted(int x[code_length], int non_sys_base_G[M][code_length - M]){
    
    for(int i = 0; i<code_length-M;i++)x[M+i] = 0;

    //int coeffs_non_sys[code_length-M];
    
    for(int i = 0;i<M;i++){
    	int coeff = rand()%Z;
    	x[i] = coeff;    
    	for (int j = 0;j<code_length-M;j++){
    		x[M+j]+=(coeff*non_sys_base_G[i][j]);
		}
	}
	
	for(int i = 0;i<code_length-M;i++)x[M+i] = x[M+i]%Z;
        
    return ;
}

/********************************************************/
//Square and multiply over Fq, already constant time
//coeff is the exponent
int square_and_multiply(int coeff){
	
	int val = 1;
	
	//find MSB
	int msb_pos = 0;
	for(int i = 0;i<LOGZ;i++){
		if(coeff&(1<<i)){
			msb_pos = i;	
		}
	}
	
	for(int i = msb_pos;i>=0;i--){
		val = (val * val)%Fq_size;
		int bit_val = ((1<<i)&coeff)>>i;
		int mult_val = 1+(E_gen-1)*bit_val;
		val = (val * mult_val)%Fq_size;	
	}
	
	return val;
}

/********************************************************/

//Compute syndrome
//Remember that e is in exponents notation
void compute_syndrome(const int e[code_length - code_dimension], const int non_sys_H[code_length - code_dimension][code_dimension], int s[code_length - code_dimension]){
		
	for(int i = 0;i<code_length-code_dimension; i++){
		int e_i = square_and_multiply(e[i]); //from exp to standard representation
		s[i] = e_i;
	}
	
	
    for(int i = 0;i<code_dimension;i++){
    	
		int e_i = square_and_multiply(e[code_length - code_dimension + i]); //from exp to standard representation
    	
    	for (int j = 0;j<code_length-code_dimension;j++){
    		s[j]=(s[j]+e_i*non_sys_H[j][i])%Fq_size;
		}
	}
	
//	for(int i = 0;i<code_length - code_dimension;i++)s[i] = s[i]%Fq_size;
      
    return ;
}

/*********************************************************************/

void compute_syndrome_Fq(const int e[code_length - code_dimension], const int non_sys_H[code_length - code_dimension][code_dimension], int s[code_length - code_dimension]){
	
	for(int i = 0;i<code_length-code_dimension; i++){
		int e_i = e[i]; //from exp to standard representation
		s[i] = e_i;
	}
	
	
    for(int i = 0;i<code_dimension;i++){
    	
		int e_i = e[code_length - code_dimension + i]; //from exp to standard representation
    	
    	for (int j = 0;j<code_length-code_dimension;j++){
    		s[j]+=(e_i*non_sys_H[j][i])%Fq_size;
		}
	}
	
	for(int i = 0;i<code_length - code_dimension;i++)s[i] = s[i]%Fq_size;
      
    return ;
}

/*********************************************************************/
//Multiply restricted objects (in exponent notation)
void multiply_restricted(const int a[code_length], const int b[code_length], int c[code_length]){
	
	for(int i = 0;i<code_length; i++){
		
		c[i] = (a[i] + b[i])%Z;
	}
	      
    return ;
}

/*********************************************************************/

//Do a * b^-1, for restricted objects (in exponent notation)
void divide_restricted(const int a[code_length], const int b[code_length], int c[code_length]){
	
	for(int i = 0;i<code_length; i++){
		
		c[i] = (a[i] - b[i]+Z)%Z;
	}
	      
    return ;
}

/*********************************************************************/

//Inverse of a restricted objects (in exponent notation)
void invert_restricted(const int a[code_length], int b[code_length]){
	
	for(int i = 0;i<code_length; i++){
		
		b[i] = Z - a[i];
	}
	      
    return ;
}

/*********************************************************************/

//Sum in Fq
void sum_Fq(const int a[code_length], const int b[code_length], int c[code_length]){
	
	
	for(int i = 0;i<code_length; i++){
		
		c[i] = (a[i] + b[i])%Fq_size;
	}
	      
    return ;
}

/*********************************************************************/

//Multiply restricted by non restricted
//a is the restricted vector (in exp notation)

void multiply_restricted_by_Fq(const int a[code_length], const int b[code_length], int c[code_length]){
	
	for(int i = 0;i<code_length; i++){
		int coeff = square_and_multiply(a[i]);
		c[i] = (coeff * b[i])%Fq_size;
	}
	      
    return ;
}


//Convert int array to string
void array2string(const int x[], const int length, char x_char[]){
	
	for(int i = 0;i<length;i++){
		x_char[2*i] = (char)((x[i] >> 8) & 0xFF);
		x_char[2*i+1] = (char)(x[i] & 0xFF);
	}
}
