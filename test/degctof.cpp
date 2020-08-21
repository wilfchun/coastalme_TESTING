#include <stdio.h>

#ifdef __cplusplus
extern"C" {
    #endif
    double DegCtoF(double [], double [], const int *);
    #ifdef __cplusplus
}
#endif


/**********************************************************************/


int main(int argc, char *argv[])
{
    const int N = 2;
    printf("C/C++ and Fortran together!\n");

    double DegreesC[N] = {32, 64};
    double DegreesF[N];

    DegCtoF(DegreesC, DegreesF, &N);
    for(int i = 0; i<N; i++){
        printf("%d : %3.1f [C] = %3.1f [F]\n", i, DegreesC[i], DegreesF[i] );
    }

    return 0;
}
