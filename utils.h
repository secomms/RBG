#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"

//Generate non systematic part of H
void gen_non_sys_H(int  non_sys_H[code_length-code_dimension][code_dimension]);
	


//Generate non systematic part for base_of_G
void gen_non_sys_base_G(int  non_sys_base_G[M ][code_length - M], int  non_sys_null_G[code_length - M][M]);

/*********************************************************/

//Verify a restricted object is in G
int verify_G(int x[code_length], const int non_sys_null_G[code_length - M][M]);

/*********************************************************************/

void gen_Fq(int x[code_length], const int seed);

/********************************************************/

//Generate restricted object
//The returned vector will be a series of exponents to assign to g
void gen_restricted(int x[code_length], int non_sys_base_G[M][code_length - M]);

/********************************************************/
//Square and multiply over Fq, already constant time
//coeff is the exponent
int square_and_multiply(int coeff);

/********************************************************/

//Compute syndrome
//Remember that e is in exponents notation
void compute_syndrome(const int e[code_length - code_dimension], const int non_sys_H[code_length - code_dimension][code_dimension], int s[code_length - code_dimension]);

/*********************************************************************/

void compute_syndrome_Fq(const int e[code_length - code_dimension], const int non_sys_H[code_length - code_dimension][code_dimension], int s[code_length - code_dimension]);

/*********************************************************************/
//Multiply restricted objects (in exponent notation)
void multiply_restricted(const int a[code_length], const int b[code_length], int c[code_length]);

/*********************************************************************/

//Do a * b^-1, for restricted objects (in exponent notation)
void divide_restricted(const int a[code_length], const int b[code_length], int c[code_length]);

/*********************************************************************/

//Inverse of a restricted objects (in exponent notation)
void invert_restricted(const int a[code_length], int b[code_length]);

/*********************************************************************/

//Sum in Fq
void sum_Fq(const int a[code_length], const int b[code_length], int c[code_length]);

/*********************************************************************/

//Multiply restricted by non restricted
//a is the restricted vector (in exp notation)

void multiply_restricted_by_Fq(const int a[code_length], const int b[code_length], int c[code_length]);


void array2string(const int x[], const int length, char x_char[]);

#endif
