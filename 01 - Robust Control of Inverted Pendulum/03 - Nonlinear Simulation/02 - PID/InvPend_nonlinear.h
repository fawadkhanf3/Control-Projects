#ifndef _InvPend_nonlinear_H_
#define _InvPend_nonlinear_H_

#define pi			(3.141592653589793)	 
#define r2d			(57.29577951308232)
#define d2r			(0.0174532925199433)
#define g0		    (9.80665)

#define Mc		    (1.0)
#define Mr			(0.25)
#define l			(0.5)
#define b			(0.05)

#define SQR(A)		((A)*(A))
#define CUB(A)		((A)*(A)*(A))
#define	MAXI(A, B)	((A) > (B) ? (A) : (B))
#define	MINI(A, B)	((A) < (B) ? (A) : (B))
#define	ssign(y)	(y < 0 ? -1 : 1)
#define U(element) 	(*uPtrs[element])

#define FLOAT_ERROR_MESSAGE(N,R) {		\
	printf("\n!%1d\t%s\n", N, R);		\
		ssSetErrorStatus(S,				\
			"DYNAM3 INIT ERROR - "		\
			R);							\
		return ;						\
	}

#define	MatMult(rowa, colb, cola, Au, Bv, Cw) \
{										\
	int i_i, j_j, k_k;					\
	double  *pCw = Cw;					\
	for(i_i=0; i_i<rowa; i_i++) {		\
		for(j_j=0; j_j<colb; j_j++) {	\
			double *pAu = Au + i_i*cola;\
			double *pBv = Bv + j_j;		\
		    *pCw = 0;					\
			for(k_k=0; k_k<cola; k_k++) { \
			  if(k_k){ ++pAu;pBv+=colb; } \
			  *pCw += *pAu * (*pBv);	\
			}							\
			++pCw;						\
		}								\
	}									\
}

#define Transpose(rowa, cola, A, B)     \
{   int k, l;                           \
    for(k=0; k<cola; k++)               \
        for(l=0; l<rowa; l++)           \
            *(*(B+l)+k) = *(*(A+k)+l);  \
}

#define MatAdd(rowa, cola, A, B, C)     \
{                                       \
    int k, l;                           \
    for(k=0; k<rowa; k++)               \
        for(l=0; l<cola; l++)           \
            *(C+k+l*rowa) = (*(A+k+l*rowa)) + (*(B+k+l*rowa)); \
}

#define MatSub(rowa, cola, A, B, C)     \
{                                       \
    int k, l;                           \
    for(k=0; k<rowa; k++)               \
        for(l=0; l<cola; l++)           \
            *(C+k+l*rowa)=(*(A+k+l*rowa))-(*(B+k+l*rowa)); \
}

#define inverse(n,A,B)		                    \
{                                               \
		int i, j, k, nr, nc;                    \
		double determinant = 1.0;               \
               nr = nc = n;                     \
      for(i=0; i < nr; ++i) {                   \
           for(j=0; j < nc; ++j) {              \
                B[i][j] = A[i][j];              \
        }                                       \
            }                                   \
		for(i=0; i < nr; ++i) {                 \
			determinant *= B[i][i];	            \
			if(determinant == 0.0) {            \
				if(B[i][i]==0.0) {              \
					printf("Inverse is not possible");           \
					break;	                                     \
				}                                                \
				else determinant = -1.0;                         \
			}                                                    \
			double r = 1.0/B[i][i];                              \
			B[i][i] = 1.0;                                       \
			for(j=0; j<nc; ++j) B[i][j] *= r;                    \
			for(k=0; k < nr; ++k) {                              \
				if(i != k && B[k][i] != 0.0) {                   \
					r = B[k][i];	                             \
					B[k][i] = 0.0;                               \
					for(j=0; j<nc; ++j)	B[k][j] -= r*B[i][j];    \
				}                                                \
			}                                                    \
		}                                                        \
	}

#endif