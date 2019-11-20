#include <iostream>

#include <cmath>
#include "isingmodel.h"
#include <iomanip>
#include <fstream>
#include <time.h>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[])
{
    double  time_start, time_end, total_time;

    // initialize MPI
    int numprocs, my_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    time_start = MPI_Wtime();
    if (argc < 8) {
        if(my_rank == 0) {
            cout << "Bad usage. Command line arguments:" << endl;
            cout << "T_initial, T_final, T_step, L, N_MCcycles, outfilename, initial state (Random or Uniform)" << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double T_initial = atof(argv[1]); // initial temperature
    double T_final   = atof(argv[2]); // final temperature
    double T_step    = atof(argv[3]); // temperature step
    int L = atoi(argv[4]);
    int MCcycles = atoi(argv[5]);
    string outfilename = argv[6];
    string initial_state = argv[7];

    // Initialize the outputfile
    ofstream outfile;
    if (my_rank == 0) {
        outfile.open(outfilename, ios::out);
        outfile << "# Initial state: " << initial_state << endl;
        outfile << "# T  -  <E>  -  <|M|>  -  C_V  -  Xsi  -  accepted" << endl;
    }
    int n_averages = 4;
    double E, M, w[17], average[n_averages], total_average[n_averages];
    // initialize lattice. First allocating memory.
    int **S;
    for (int i=0; i<L; i++) {
        S = new int *[L];
    }
    for (int i=0; i<L; i++) {
        S[i] = new int [L];
    }


    /*
    Determine number of intervall which are used by all processes
    myloop_begin gives the starting point on process my_rank
    myloop_end gives the end point for summation on process my_rank
    */

    int no_of_intervals = MCcycles/numprocs;
    int myloop_begin = my_rank*no_of_intervals + 1;
    int myloop_end = (my_rank + 1)*no_of_intervals;
    if ( (my_rank == numprocs-1) && ( myloop_end < MCcycles) ) myloop_end = MCcycles;

    // Now let all nodes know about the common variables
    // They are: L, T_initial, T_final, T_step
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_initial, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&T_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    long idum = -1-my_rank; // starting point for the RNG

    initialize(S, E, M, L, idum, initial_state);
    for (double T = T_initial; T <= T_final; T+=T_step) {
        // initializing the lattice with random or uniform spin configuration

        // pre-calculating the possible changes in energy

        for (int deltaE = -8; deltaE <= 8; deltaE++) {
            w[deltaE+8] = 0;
        }
        for (int deltaE = -8; deltaE <= 8; deltaE+=4) {
            w[deltaE+8] = exp(-deltaE/T);
        }

        for (int i = 0; i < n_averages; i++) average[i] = 0.0;
        for (int i = 0; i < n_averages; i++) total_average[i] = 0.0;

        int accepted = 0; // number of accepted energy states
        for (int cycles = myloop_begin; cycles <= myloop_end; cycles++) {
            Metropolis(S, E, M, w, L, idum, accepted);
            // update expectation values
            average[0] += E;    average[1] += E*E;
            average[2] += fabs(M); average[3] += fabs(M)*fabs(M);
        }

        // find total average
        for (int i = 0; i < n_averages; i++) {
            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }


        double norm = 1/((double) (MCcycles));  // divided by total number of cycles
        double E_average = total_average[0]*norm;
        double E2_average = total_average[1]*norm;
        double M_abs_average = total_average[2]*norm;
        double M_abs2_average = total_average[3]*norm;
        // all expectation values are per spin, divide by L*L
        double E_variance = (E2_average- E_average*E_average)/L/L;
        double M_abs_variance = (M_abs2_average - M_abs_average*M_abs_average)/L/L;

        outfile << setw(15) << setprecision(8) << T;
        outfile << setw(15) << setprecision(8) << E_average/L/L;
        outfile << setw(15) << setprecision(8) << M_abs_average/L/L;
        outfile << setw(15) << setprecision(8) << E_variance/T/T; // Divide by T*T because heat capacity
        outfile << setw(15) << setprecision(8) << M_abs_variance/T; // Divide by T because susceptibility
        outfile << setw(15) << setprecision(8) << accepted/((double) MCcycles) << endl;

        if (my_rank == 0) cout << "Finished " << T << " of " << T_final << endl;
    }
    time_end = MPI_Wtime();
    total_time = time_end-time_start;

    if (my_rank == 0) {
        cout << "Time: " << total_time << " on " << numprocs << " processors" << endl;
    }

    delete[] S;
    outfile.close();
    MPI_Finalize();
    return 0;
}
