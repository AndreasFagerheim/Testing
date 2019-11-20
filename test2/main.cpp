//task d)



#include <iostream>
#include "../lib.h"
#include <cmath>
#include "isingmodel.h"
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 6) {
        cout << "Bad usage. Command line arguments:" << endl;
        cout << "T, L, MCcylces, outfilename, initial state (Random or Uniform)" << endl << endl;
    }

    double T = atof(argv[1]); // Temperature
    int L = atoi(argv[2]);
    int MCcycles = atoi(argv[3]);
    string outfilename = argv[4];
    string initial_state = argv[5];

    long idum = -1; // starting point for the RNG
    double E, M;

    // pre-calculating the possible changes in energy
    double w[17];
    double average[5];

    for (int deltaE = -8; deltaE <= 8; deltaE++) {
        w[deltaE+8] = 0;
    }
    for (int deltaE = -8; deltaE <= 8; deltaE+=4) {
        w[deltaE+8] = exp(-deltaE/T);
    }

    // initialize lattice
    int **S;
    for (int i=0; i<L; i++) {
        S = new int *[L];
    }
    for (int i=0; i<L; i++) {
        S[i] = new int [L];
    }
    // initializing the lattice with random or uniform spin configuration
    initialize(S, E, M, L, idum, initial_state);


    ofstream outfile;
    outfile.open(outfilename, ios::out);
    outfile << "# Initial state: " << initial_state << endl;
    outfile << "# Temperature: " << T << endl;
    outfile << "# Lattice size L = " << L << endl;
    outfile << "# MCcycles  -  <E>/particle  -  <|M|>/particle  -  E_total  -  N_accepted" << endl;


    for (int i = 0; i < 5; i++) average[i] = 0.0;

    E = calc_total_energy(S, L);
    double norm, Eaverage, Mabsaverage;
    int accepted = 0; // number of accepted energy states
    for (int cycles=1; cycles<=MCcycles; cycles++) {
        Metropolis(S, E, M, w, L, idum, accepted);
        // update expectation values
        average[0] += E;    average[1] += E*E;
        average[2] += M;    average[3] += M*M;  average[4] += fabs(M);

        norm = 1/((double) cycles);
        Eaverage = average[0]*norm;
        Mabsaverage = average[4]*norm;

        outfile << setw(15) << setprecision(8) << cycles;
        outfile << setw(15) << setprecision(8) << Eaverage/L/L;     // <E>
        outfile << setw(15) << setprecision(8) << Mabsaverage/L/L;  // <|M|>
        outfile << setw(15) << setprecision(8) << calc_total_energy(S, L);
        outfile << setw(15) << setprecision(8) << accepted/((double) cycles) << endl;
    }
    double Evariance = (average[1]/MCcycles - average[0]*average[0]/MCcycles/MCcycles)/L/L;
    cout << "Variance = " << Evariance << endl;
    outfile.close();

    delete[] S;
    return 0;
}
