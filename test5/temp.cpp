/* Program for simulating a 2D lattice
 * with the Ising model by using the
 * Metropolis algorithm
 *
 * THIS PROGRAM IS CUSTOMIZED TO RUN ON LAB MACHINES
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "lib.h"
#include "time.h"
#include "mpi.h"

using namespace std;

/* 1. Create a state, compute the difference
 *      - a matrix rep. one state
 * 2. Flip one spin, compute the energy
 *      - choose a random site in the lattice by using Monte Carlo
 * 3. Compute the energy difference
 * 4. If E_delta <= 0, go to 7.
 * 5. Calc. w = exp(-beta*E_delta)
 * 6. Compare w and random number r, if r<=w, keep the new state
 * 7. Update various expectation values
 * 8. Repeat 2-7 a given number of times
 *
 * b)   initial_temperature = 1.0
 *      L = 2 (i.e. four spins)
 */

// ensuring boundary conditions by
// "coupling togheter" the elements
// on the edges
inline int periodic(int i, int dimension, int add) {
    return (i+dimension+add) % (dimension);
    // The add-variable is either 1 or -1 depending on which side
    // of the element we are on, e.g.
    // state_matrix[MC_row][periodic(MC_column, L_spins, -1)] is
    // the element below the one we're looking at
}

void metropolis_algo(int **state_matrix, int L_spins, double &energy, double &magnetization, double *w, long &idum, int &accepted, int **neighbour);
void CreateStartingState(int **, double, double , double , long &, bool);
void ProbabilityForGivenEnergy(int **state_matrix, int L_spins, double &energy, int **neighbour);
void SavingResultsForDifferentTemperatures(string FileName, double *e_var, double *m_var, double *mean_mag, double *temperatures, int no_of_steps, double *avg_energy);
void ExpectationValues(double *values, double normalization, double &accepted_total, double &HeatCapacity, double &Susceptibility, double &MeanMagnetization, int L, double &avg_energy);


int main(int argc, char* argv[]){
    //------------------------------------------------------
    // Read in initial parameters on command line
    double start_temp = atof(argv[1]);// 2.1; // argv[1]
    double final_temp = atof(argv[2]);//2.6; // argv[2]
    double temp_step  = atof(argv[3]);//0.01; // argv[3]
    int L = atoi(argv[4]); //80; // argv[4] // Size of lattice 40,60,80
    double MC = atof(argv[5]);//1000000; //argv[5]
    string filename = argv[6];

    // Try to give values which give a N which is a multiple of 8 (or 4)
    int N = (final_temp - start_temp)/temp_step + 1;
    bool random = true; // Starting with a random or non-random state

    //------------------------------------------------------
    //  MPI initializations
    int my_rank, numprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Assigning start and stop temperatures for each thread
    int no_intervalls = N/numprocs;
    int myloop_begin = my_rank*no_intervalls;
    int myloop_end = (my_rank+1)*no_intervalls -1;
    if ( (my_rank == numprocs-1) && ( myloop_end <= N) ) myloop_end = N-1;
    // broadcast to all nodes common variables
    MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&MC, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //------------------------------------------------------
    double values[5]; // Array for keeping our results
    double energy=0, magnetization=0;
    double w[17]; // Array containing all possible energies for the system kind of sth what

    // Creating array which contain information about the neighbouring
    // spins when we use periodic boundary conditions.
    int **neighbour;
    neighbour = (int **) matrix(L,2,sizeof(L));
    for(int i=0;i<L;i++){
        neighbour[i][0] = periodic(i,L,-1);
        neighbour[i][1] = periodic(i,L,1);
    }

    // Arrays for keeping results from each thread
    double HeatCapacity[N];
    double Susceptibility[N];
    double MeanMagnetization[N];
    double temperature[N];
    double accepted_total[N];
    double avg_energy[N];
    for(int i=0;i<N;i++) temperature[i] = start_temp + i*temp_step;

    // Allocating memory for state matrix and then filling it
    int**state_matrix;

    // Arrays for containing merged results from the threads
    double tot_HeatCapacity[N];
    double tot_Susceptibility[N];
    double tot_MeanMagnetization[N];
    double tot_temperature[N];
    double tot_avg_energy[N];

    // Setting array elements to zero to avoid fuck ups
    for(int q=0;q<N;q++) tot_HeatCapacity[q]=0,tot_Susceptibility[q]=0,tot_MeanMagnetization[q]=0,tot_temperature[q]=0, tot_avg_energy[q]=0;
    for(int q=0;q<N;q++) HeatCapacity[q]=0,Susceptibility[q]=0,MeanMagnetization[q]=0, accepted_total[q]=0, avg_energy[q]=0;

    // Ensuring that each thread get different seed
    long idum = -1 - my_rank;

    // Each thread has its own state matrix
    state_matrix = (int **) matrix(L, L, sizeof(L));
    CreateStartingState(state_matrix, L, energy=0, magnetization=0, idum, random);

    int q = 0;
    // Running one temperature on each thread in order to get as many MC-samples
    // of each temperature.
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        q = cycles;
        // Setting values to zero
        for(int i=0;i<5;i++) values[i] = 0;
        // Precalculating the different temperature differences
        for( int delta_E =-8; delta_E <= 8; delta_E++) w[delta_E+8] = 0;
        for( int delta_E =-8; delta_E <= 8; delta_E+=4) w[delta_E+8] = exp(-delta_E/temperature[q]);
        for(int n=0;n<MC; n++){
            int accepted = 0;
            metropolis_algo(state_matrix, L, energy, magnetization, w, idum, accepted, neighbour);
            values[0] += energy;
            values[1] += energy*energy;
            values[2] += magnetization;
            values[3] += magnetization*magnetization;
            values[4] += fabs(magnetization);
            accepted_total[q] += accepted;
        }
        ExpectationValues(values, MC, accepted_total[q], HeatCapacity[q], Susceptibility[q], MeanMagnetization[q], L, avg_energy[q]);
        if (my_rank == 0){
            cout << "Nå har jeg kommet til: " << temperature[q] << endl;
        }
    }
    free_matrix(((void **) state_matrix)); // free memory

    // Merging results from every thread
    MPI_Reduce(&HeatCapacity, &tot_HeatCapacity, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Susceptibility, &tot_Susceptibility, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&MeanMagnetization, &tot_MeanMagnetization, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&temperature, &tot_temperature, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_energy, &tot_avg_energy, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // End MPI
    MPI_Finalize ();

    // Only one thread should write the file!
    if ( my_rank == 0) {
        for(int q=0;q<N;q++) cout << tot_HeatCapacity[q] << " ";
        SavingResultsForDifferentTemperatures(filename, tot_HeatCapacity, tot_Susceptibility, tot_MeanMagnetization, temperature, N, tot_avg_energy);
    }

    return 0;
}

void metropolis_algo(int **state_matrix, int L_spins, double &energy, double &magnetization, double *w, long &idum, int &accepted, int **neighbour){
    for(int row=0;row<L_spins;row++){
        for(int column=0;column<L_spins;column++){
            // Picking random site in lattice
            int MC_row =    (int) (ran2(&idum)*(double)L_spins);
            int MC_column = (int) (ran2(&idum)*(double)L_spins);
            // calculating energy difference if one spin is flipped
            int delta_E = 2*state_matrix[MC_row][MC_column]*
                    (state_matrix[MC_row][neighbour[MC_column][0]] +
                    state_matrix[neighbour[MC_row][0]][MC_column] +
                    state_matrix[MC_row][neighbour[MC_column][1]] +
                    state_matrix[neighbour[MC_row][1]][MC_column]);
            // comparing random number with energy difference
            if( delta_E <= 0 or ran2(&idum) <= w[delta_E + 8]){
                accepted += 1;
                state_matrix[MC_row][MC_column] *= -1;
                energy += (double) delta_E;
                magnetization += (double) 2*state_matrix[MC_row][MC_column];
            }
        }
    }
}


void CreateStartingState(int **state_matrix, double L_spins, double energy, double magnetization, long &idum, bool random){
    if(random){
        // Filling state matrix with random spins
        for(int row=0;row<L_spins;row++){
            for(int column=0;column<L_spins;column++){
                if(ran2(&idum) <= 0.5) state_matrix[row][column] = 1.;
                else state_matrix[row][column] = -1;
            }
        }
    }

    else{
        // Filling state matrix with spin one states
        for(int row=0;row<L_spins;row++){
            for(int column=0;column<L_spins;column++){
                state_matrix[row][column] = 1.;

            }
        }
    }
    // Calculating energy and magnetization of initial state
    for(int row=0;row<L_spins;row++){
        for(int column=0;column<L_spins;column++){
            energy -=  (double) state_matrix[row][column]*(state_matrix[periodic(row,L_spins,-1)][column] +
                    state_matrix[row][periodic(column, L_spins, -1)]);
            magnetization += (double) state_matrix[row][column];
        }
    }

}

void ProbabilityForGivenEnergy(int **state_matrix, int L_spins, double &energy, int **neighbour){
    // Function for calculating the probability of a given energy
    energy = 0;
    for(int row=0;row<L_spins;row++){
        for(int column=0;column<L_spins;column++){
            energy -= (double) state_matrix[row][column]*(state_matrix[neighbour[row][0]][column] +
                    state_matrix[row][neighbour[column][0]]);
        }
    }
}

void SavingResultsForDifferentTemperatures(string FileName, double *e_var, double *m_var, double *mean_mag,
                                           double *temperatures, int no_of_steps, double *avg_energy){
    ofstream myfile;
    myfile.open (FileName);
    myfile << "e_var " << "= [";
    for (int i=0; i<no_of_steps; i++){
        myfile << e_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "m_var " << "= [";
    for (int i=0; i<no_of_steps; i++){
        myfile << m_var[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "mean_mag " << "= [";
    for (int i=0; i<no_of_steps; i++){
        myfile << mean_mag[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "temperatures " << "= [";
    for (int i=0; i<no_of_steps; i++){
        myfile << temperatures[i] << ", ";
    }
    myfile << "];" << endl;

    myfile << "avg_energy " << "= [";
    for (int i=0; i<no_of_steps; i++){
        myfile << avg_energy[i] << ", ";
    }
    myfile << "];" << endl;

    myfile.close();
}

void ExpectationValues(double *values, double normalization, double &accepted_total, double &HeatCapacity, double &Susceptibility, double &MeanMagnetization, int L, double &avg_energy){
    accepted_total = accepted_total/normalization;
    HeatCapacity = (values[1]/normalization - values[0]*values[0]/(normalization*normalization)) / (L*L);
    Susceptibility = (values[3]/normalization - values[4]*values[4]/(normalization*normalization))/(L*L);
    MeanMagnetization = values[4]/(normalization*L*L);
    avg_energy = values[0]/(normalization*L*L);
}
