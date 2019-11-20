

//chekk if uniform or randim init state.
#include "isingmodel.h"
#include <iostream>
#include <cmath>

using namespace std;

double ran2(long *idum);
int periodic(int index, int size, int add)
{
    return (index + size + add) % size;
}


void Metropolis(int **S, double &E, double &M, double *w, int L, long &idum, int &accepted)
{
    // loop over all the elements in the lattice
    for (int y=0; y < L; y++) {
        for (int x=0; x < L; x++) {
            // Choose two random numbers in [0, L-1]. This is where I position myself in the lattice
            int ix = (int) (ran2(&idum)*(double)L);
            int iy = (int) (ran2(&idum)*(double)L);
            // calculate the change in energy of this state
            int deltaE = 2*S[iy][ix]*
                    (S[iy][periodic(ix, L, 1)] +
                    S[iy][periodic(ix, L, -1)] +
                    S[periodic(iy, L, 1)][ix] +
                    S[periodic(iy, L, -1)][ix]);
            // if deltaE <= 0 accept the new configuration, if not, check if a random number <= exp(deltaE/T) and accept
            if (ran2(&idum) <= w[deltaE+8] ) {
                S[iy][ix] *= -1;
                accepted += 1;
                E += (double) deltaE;
                M += (double) 2*S[iy][ix];
            }
        }
    }
}



double calc_total_energy(int **S, int L)
{
    double E = 0;
    for (int y = 0; y < L; y++) {
        for (int x = 0; x < L; x++) {
            E -= (double) S[y][x]*
                    (S[periodic(y,L,-1)][x] +
                    S[y][periodic(x, L, -1)]);
        }
    }
    return E;
}

void initialize(int **S, double &E, double &M, int L, long &idum, string initial_state)
{


    E = M = 0.0;
    if (initial_state == "Random") {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                S[i][j] = round(ran2(&idum))*2 - 1;
                M += (double) S[i][j];
            }
        }
    }
    if (initial_state == "Uniform") {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                S[i][j] = 1;
                M += (double) S[i][j];
            }
        }
    }

    E = calc_total_energy(S, L);

}
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
   int            j;
   long           k;
   static long    idum2 = 123456789;
   static long    iy=0;
   static long    iv[NTAB];
   double         temp;

   if(*idum <= 0) {
      if(-(*idum) < 1) *idum = 1;
      else             *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ1;
     *idum = IA1*(*idum - k*IQ1) - k*IR1;
     if(*idum < 0) *idum +=  IM1;
     if(j < NTAB)  iv[j]  = *idum;
      }
      iy=iv[0];
   }
   k     = (*idum)/IQ1;
   *idum = IA1*(*idum - k*IQ1) - k*IR1;
   if(*idum < 0) *idum += IM1;
   k     = idum2/IQ2;
   idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
   if(idum2 < 0) idum2 += IM2;
   j     = iy/NDIV;
   iy    = iv[j] - idum2;
   iv[j] = *idum;
   if(iy < 1) iy += IMM1;
   if((temp = AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
