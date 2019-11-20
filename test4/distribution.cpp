#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <boost/timer.hpp>
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
void InitializeLattice_RNG(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &);
void MetropolisSampling_Terminate(int, int, double);
// prints to file the results of the calculations
void WriteResultstoFile(int, int, double, vec);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename;
  int NSpins, MCcycles;
  double InitialTemp, FinalTemp, TempStep;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
    exit(1);
  }
  if (argc > 1) {
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MCcycles = atoi(argv[3]);
    InitialTemp = atof(argv[4]);
    FinalTemp = atof(argv[5]);
    TempStep = atof(argv[6]);
  }
  // Declare new file name and add lattice size to file name
  string fileout = filename;
  string argument = to_string(NSpins);
  fileout.append("N"+argument);
  fileout.append(".txt");
  ofile.open(fileout);
  ofile << "T     E/N       ";
  ofile << "E-var/TT       "<< endl;
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
    MetropolisSampling_Terminate(NSpins, 1000000, Temperature);
    vec ExpectationValues = zeros<mat>(5);
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues);
    //
  }
  ofile.close();  // close output file
  return 0;
}



// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, vec &ExpectationValues)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values



  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17);



  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
      vec values = zeros<mat>(2);
      values(0) = Energy*MCcycles;
      values(1) = Energy*Energy*MCcycles;
      WriteResultstoFile(NSpins, MCcycles, Temperature, values);
      for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){
          int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
          int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
          int deltaE =  2*SpinMatrix(ix,iy)*
            (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
             SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
             SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
             SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
          if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
            SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
            MagneticMoment += (double) 2*SpinMatrix(ix,iy);
            Energy += (double) deltaE;
          }


        }
      }


    // update expectation values  for local node
  }
} // end of Metropolis sampling over spins


void MetropolisSampling_Terminate(int NSpins, int MCcycles, double Temperature)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values



  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17);



  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles++){
      // The sweep over the lattice, looping over all spin sites
      for(int x =0; x < NSpins; x++) {
        for (int y= 0; y < NSpins; y++){


          int ix = (int) (RandomNumberGenerator(gen)*(double)NSpins);
          int iy = (int) (RandomNumberGenerator(gen)*(double)NSpins);
          int deltaE =  2*SpinMatrix(ix,iy)*
            (SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
             SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
             SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
             SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
          if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
            SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
            MagneticMoment += (double) 2*SpinMatrix(ix,iy);
            Energy += (double) deltaE;
          }


        }
      }
  }
} // end of Metropolis sampling over spins



// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  // setup spin matrix and initial magnetization
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }
  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
    (SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
     SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise

// function to initialise energy, spin matrix and magnetization
void InitializeLattice_RNG(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  // setup spin matrix and initial magnetization
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      if (RandomNumberGenerator(gen) > 0.5){
        SpinMatrix(x,y) = 1.0;
      } else {
        SpinMatrix(x,y) = -1.0;
      }
      MagneticMoment +=  (double) SpinMatrix(x,y);
    }
  }
  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
    (SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
     SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialise


void WriteResultstoFile(int NSpins, int MCcycles, double temperature, vec ExpectationValues)
{
  double norm = 1.0/((double) (MCcycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/NSpins/NSpins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/NSpins/NSpins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/NSpins/NSpins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues;
  ofile << setw(15) << setprecision(8) << E2_ExpectationValues << endl;
} // end output function
