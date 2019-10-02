
#include "constants.h"
#include "functions.h"


int main() {
    //Initialize variables
    double theta[L][L][L], phi[L][L][L], E, T = T_init, E_sum, E_sumsquares;
    int near_n[L][2];


    //store locations of nearest neighbours using PBCs
    get_neighbours(near_n);

    //Initialize the spin directions for each lattice index
    //random_state(theta, phi);
    uniform_state(theta, phi);


    //calculate initial net magnetisation and output to file
    magnetization(theta, phi);

    //Calculate total system energy and output to file
    total_energy(theta, phi, near_n, E);


    //MC loops inside a loop that varies temperature
    for(int b = 0; b < T_intervals; b++) {
        MC_steps(theta, phi, near_n, T, E, E_sum, E_sumsquares);

        //Calculates and outputs the heat capacity
        heat_cap(E_sum, E_sumsquares, T);

        magnetization(theta, phi);

        total_energy(theta, phi, near_n, E);

        //Outputs T to file then adds T_step
        iterate_T(T);
    }
//end
}
