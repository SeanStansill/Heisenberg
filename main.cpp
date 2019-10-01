#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include "constants.h"
#include "functions.h"


int main() {
    //Initialize variables
    double theta[L][L][L], phi[L][L][L], theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {}, E, T = T_init, C, E1, E2, dE, sum, sq_sum;
    int near_n[L][2];
    std::vector<double> M;

    //Seed rng
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> uniform_distribution;
    std::uniform_real_distribution<double> decimal(0.0, 1.0);
    std::uniform_real_distribution<double> r_angle(0.0, 2.0*pi);


    //store locations of nearest neighbours using PBCs
    get_neighbours(near_n);

    //Initialize the spin directions for each lattice index
    //random_state(theta, phi);
    uniform_state(theta, phi);


    //calculate initial net magnetisation
    M = magnetization(theta, phi);

    write_M(M);
    E = total_energy(theta, phi, near_n);
    write_E(E);


    int i, j, k;

    //MC loops inside a loop that varies temperature
    for(int b = 0; b < T_intervals; b++) {
        sum = 0.0;
        sq_sum = 0.0;
        for (int a = 0; a < nsteps; a++) {
            i = uniform_distribution(mt) % L;
            j = uniform_distribution(mt) % L;
            k = uniform_distribution(mt) % L;

            E1 = local_energy(i, j, k, theta, phi, near_n);

            theta_trial[i][j][k] = 2*pi*decimal(mt);
            phi_trial[i][j][k] = acos((2*decimal(mt))-1);

            E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

            dE = E2 - E1;
            if(decimal(mt) <= exp(std::min(0.0, (((-J)/(k_B*T)) * dE)))){
                theta[i][j][k] = theta_trial[i][j][k];
                phi[i][j][k] = phi_trial[i][j][k];
            }

            E = total_energy(theta, phi, near_n);
            sum += E;
            sq_sum += (E*E);

        }

        //calculate and output the magnetisation after every set of MC loops
        sum = sum/nsteps;
        sq_sum = sq_sum/nsteps;
        C = heat_cap(sum, sq_sum, T);
        write_C(C);

        M = magnetization(theta, phi);

        write_M(M);

        write_E(E);

        write_T(T);

        T += T_step;
    }
//end
}
