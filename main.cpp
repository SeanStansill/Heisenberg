#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include "constants.h"
#include "functions.h"


int main() {
    //Initialize variables
    double theta[N][N][N], phi[N][N][N], theta_trial[N][N][N] = {}, phi_trial[N][N][N] = {}, E, T = T_init, e_squared, C, E1, E2, dE;
    int near_n[N][2];
    std::vector<double> M;

    //Seed rng
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> uniform_distribution;;
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


    int i, j, k;
    for (int a = 0; a < nsteps; a++) {
        i = uniform_distribution(mt) % N;
        j = uniform_distribution(mt) % N;
        k = uniform_distribution(mt) % N;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = r_angle(mt);
        phi_trial[i][j][k] = r_angle(mt);

        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if (decimal(mt) <= exp(std::min(0.0, ((-J / (k_B * T)) * dE)))) {
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }
    }

    //MC loops inside a loop that varies temperature
    for(int b = 0; b < T_intervals; b++) {
        for (int a = 0; a < nsteps; a++) {
            i = uniform_distribution(mt) % N;
            j = uniform_distribution(mt) % N;
            k = uniform_distribution(mt) % N;

            E1 = local_energy(i, j, k, theta, phi, near_n);

            theta_trial[i][j][k] = 2*pi*decimal(mt);
            phi_trial[i][j][k] = acos((2*decimal(mt))-1);

            E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

            dE = E2 - E1;
            if(decimal(mt) <= exp(std::min(0.0, ((-J/(k_B*T)) * dE)))){
                theta[i][j][k] = theta_trial[i][j][k];
                phi[i][j][k] = phi_trial[i][j][k];
            }
        }

        //calculate and output the magnetisation after every set of MC loops
        M = magnetization(theta, phi);
        E = total_energy(theta, phi, near_n);
        //e_squared = sum_squares(theta, phi, near_n);
        write_E(E);

        C = heat_cap(theta, phi, near_n, T);

        write_C(C);

        write_M(M);

        write_T(T);

        T += T_step;
    }
//end
}