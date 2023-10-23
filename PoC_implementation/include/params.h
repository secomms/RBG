#define Fq_size              971    //finite field size
#define E_gen                  4    //generator of multiplicative subgroup
#define Z                     97    //size of mult subgroup
#define CODE_LENGTH_N         44    //code length
#define CODE_DIMENSION_K      26    //code dimension
#define M                     26    //dimension of basis for G
#define N                     32   //number of reps per round
#define T                     42    //number of rounds


#define SEC_LEVEL           128   
#define SEEDBYTES           (SEC_LEVEL/8)    //Number of bytes for seed
#define HASHBYTES           (SEEDBYTES*2)    //Number of bytes in hash
#define NUMBYTES            114    //Number of bytes for sigma + v

#define IS_REPRESENTABLE_IN_D_BITS(D, N)                \
  (((unsigned long) N>=(1UL << (D-1)) && (unsigned long) N<(1UL << D)) ? D : -1)
#define BITS_TO_REPRESENT(N)                            \
  (N == 0 ? 1 : (15                                     \
                 + IS_REPRESENTABLE_IN_D_BITS( 1, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 2, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 3, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 4, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 5, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 6, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 7, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 8, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS( 9, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(10, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(11, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(12, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(13, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(14, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(15, N)    \
                 + IS_REPRESENTABLE_IN_D_BITS(16, N)    \
                 )                                      \
   )

#define LOGZ          BITS_TO_REPRESENT(Z)    //Number of bits to represent Z
#define LOGQ          BITS_TO_REPRESENT(Q)
#define LOGN          ( (BITS_TO_REPRESENT(N) > BITS_TO_REPRESENT(N-1)) ? BITS_TO_REPRESENT(N-1) : BITS_TO_REPRESENT(N) )
  
//Comment if you don't want to print values (useful for debugging)
/* #define PRINTING               1 */
