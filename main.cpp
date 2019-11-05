
#include "constants.h"
#include "functions.h"
#include <vector>
#include <cmath>

int main() {

    double theta[L][L][L], phi[L][L][L], E, T = T_init, E_sum, E_sumsquares, C, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, chi, dSx, dSy, dSz, dS;
    int near_n[L][2];

    write_N(N);

    get_nearest_neighbours(near_n);

    //random_spin_state(theta, phi);
    uniform_spin_state(theta, phi);


    magnetization(theta, phi, Mx, My, Mz, M);
    write_M(Mx, My, Mz, M);

    total_energy(theta, phi, near_n, E);
    write_E(E);
    //thermalize_typewriter(theta, phi, near_n, T); // For random state, initial thermalisation will take much longer than after reaching saturation
    //Note to Old Self: It doesn't take that long as 50000 sweeps is sufficient to ensure the state is saturated at low T

    for (int b = 0; b < T_intervals; b++) {
        thermalize_typewriter(theta, phi, near_n, T);

        //MC_loop(theta, phi, near_n, T, E, E_sum, E_sumsquares, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, b, dSx, dSy, dSz, dS);
        MC_typewriter(theta, phi, near_n, T, E, E_sum, E_sumsquares, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, b, dSx, dSy, dSz, dS);

        heat_cap(C, E_sum, E_sumsquares, T);
        write_C(C);
        susceptibility(chi, M_sum, M_sumsquares, T);
        write_chi(chi);

        magnetization(theta, phi, Mx, My, Mz, M);
        write_M(Mx, My, Mz, M);

        write_E(E);

        iterate_T(T);
    }
}

//fftw instead of autocorrelation