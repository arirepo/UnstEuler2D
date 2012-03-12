#ifndef _INTEGRAL_OP_H_
#define _INTEGRAL_OP_H_


int S_op(int colength,int z, int *NBSCHTYPE, int len_NBSCHTYPE, int INTSCHTYPE, double **S, double **f1_mod, short filter);

int S_safe(int colength,int z, int *NBSCHTYPE, int len_NBSCHTYPE, int INTSCHTYPE, double **S, double **f1_mod);

void test_S(void);

#endif
