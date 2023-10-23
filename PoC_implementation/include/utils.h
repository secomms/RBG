#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "params.h"

//Generate non systematic part of H
void gen_non_sys_H(int
                   non_sys_H[CODE_LENGTH_N - CODE_DIMENSION_K][CODE_DIMENSION_K]);



//Generate non systematic part for base_of_G
void gen_non_sys_base_G(int non_sys_base_G[M][CODE_LENGTH_N - M],
                        int non_sys_null_G[CODE_LENGTH_N - M][M]);

/*********************************************************/

//Verify a restricted object is in G
int verify_G(int x[CODE_LENGTH_N],
             const int non_sys_null_G[CODE_LENGTH_N - M][M]);

/*********************************************************************/

void gen_Fq(int x[CODE_LENGTH_N], const int seed);

/********************************************************/

//Generate restricted object
//The returned vector will be a series of exponents to assign to g
void gen_restricted(int x[CODE_LENGTH_N],
                    int non_sys_base_G[M][CODE_LENGTH_N - M]);

/********************************************************/
//Square and multiply over Fq, already constant time
//coeff is the exponent
int square_and_multiply(int coeff);

/********************************************************/

//Compute syndrome
//Remember that e is in exponents notation
/* void compute_syndrome(const int e[CODE_LENGTH_N - CODE_DIMENSION_K], */
void compute_syndrome(const int e[CODE_LENGTH_N],
                      const int non_sys_H[CODE_LENGTH_N -
                            CODE_DIMENSION_K][CODE_DIMENSION_K],
                      int s[CODE_LENGTH_N - CODE_DIMENSION_K]);

/*********************************************************************/

/* void compute_syndrome_Fq(const int e[CODE_LENGTH_N - CODE_DIMENSION_K], */
void compute_syndrome_Fq(const int e[CODE_LENGTH_N],
                         const int non_sys_H[CODE_LENGTH_N -
                                           CODE_DIMENSION_K]
                         [CODE_DIMENSION_K],
                         int s[CODE_LENGTH_N - CODE_DIMENSION_K]);

/*********************************************************************/
//Multiply restricted objects (in exponent notation)
void multiply_restricted(const int a[CODE_LENGTH_N],
                         const int b[CODE_LENGTH_N], int c[CODE_LENGTH_N]);

/*********************************************************************/

//Do a * b^-1, for restricted objects (in exponent notation)
void divide_restricted(const int a[CODE_LENGTH_N], const int b[CODE_LENGTH_N],
                       int c[CODE_LENGTH_N]);

/*********************************************************************/

//Inverse of a restricted objects (in exponent notation)
void invert_restricted(const int a[CODE_LENGTH_N], int b[CODE_LENGTH_N]);

/*********************************************************************/

//Sum in Fq
void sum_Fq(const int a[CODE_LENGTH_N], const int b[CODE_LENGTH_N],
            int c[CODE_LENGTH_N]);

/*********************************************************************/

//Multiply restricted by non restricted
//a is the restricted vector (in exp notation)

void multiply_restricted_by_Fq(const int a[CODE_LENGTH_N],
                               const int b[CODE_LENGTH_N],
                               int c[CODE_LENGTH_N]);


void array2string(const int x[], const int length, unsigned char x_char[]);

#endif
