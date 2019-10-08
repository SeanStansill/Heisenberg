
#include "constants.h"
#include "functions.h"
#include <vector>
#include <cmath>

int main() {
    double E_initial = 0;
    double zero_timeE[steps] = {};
    double timeE[steps] = {};
    double ExpectationZero_T;
    double Corr = 10;

    double theta[L][L][L], phi[L][L][L], E, T = T_init, E_sum, E_sumsquares, C, Mx, My, Mz, sigma;
    int near_n[L][2];

    get_nearest_neighbours(near_n);

    //random_spin_state(theta, phi);
    uniform_spin_state(theta, phi);


    magnetization(theta, phi, Mx, My, Mz);
    write_M(Mx, My, Mz);

    total_energy(theta, phi, near_n, E);
    write_E(E);


    for (int b = 0; b < T_intervals; b++) {
        //thermalize(theta, phi, near_n, T);

        //initial_corr(theta, phi, near_n, T, E, E_sum, E_sumsquares, E_initial, zero_timeE, timeE, ExpectationZero_T);
        //MC_loopautocorr(theta, phi, near_n, T, E, E_sum, E_sumsquares, E_initial, zero_timeE, timeE, ExpectationZero_T, Corr);

        MC_loop(theta, phi, near_n, T, E, E_sum, E_sumsquares);

        heat_cap(C, E_sum, E_sumsquares, T);
        write_C(C);

        magnetization(theta, phi, Mx, My, Mz);
        write_M(Mx, My, Mz);

        total_energy(theta, phi, near_n, E);
        write_E(E);

        std::cout << b << " " << sqrt((Mx*Mx) + (My*My) + (Mz*Mz)) << std::endl;

        iterate_T(T);
    }
}