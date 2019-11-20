#include <iostream>
using namespace std;

int periodic(int index, int size, int add);



void Metropolis(int** S, double &E, double &M, double *w, int L, long &idum, int &accepted);

double calc_total_energy(int **S, int L);

void initialize(int **S, double &E, double &M, int L, long &idum, string initial_state);
