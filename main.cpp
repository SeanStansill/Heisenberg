
#include "constants.h"
#include "functions.h"
#include <vector>
#include <cmath>
#include <omp.h>
#include <random>
#include <random>
#include <iostream>

int main() {


    double theta[N], phi[N], E, T = T_init, E_sum, E_sumsquares, C, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, chi, moves_accepted;
    int up[N], down[N], left[N], right[N], backwards[N], forwards[N];

    write_N(N);

    get_nearest_neighbours(up, down, left, right, backwards, forwards);

    //random_spin_state(theta, phi);


    int n = omp_get_num_procs();
    omp_set_num_threads(n/2); // n-1 so that when running locally the OS is still snappy and can access a decent amount of power
#pragma omp parallel default(none) shared(up, down, left, right, backwards, forwards) private(theta, phi, E, E_sum, E_sumsquares, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, chi, T, C, moves_accepted)
{
#pragma omp for schedule(dynamic)
        for (int b = 0; b < T_intervals; b++) {
            random_spin_state(theta, phi);
            T = T_init + (b * T_step);
            thermalize_parallel(theta, phi, up, down, left, right, backwards, forwards, T);

            magnetization(theta, phi, Mx, My, Mz, M);
            total_energy(theta, phi, up, down, left, right, backwards, forwards, E);

            MC_parallel(theta, phi, up, down, left, right, backwards, forwards, T, E, E_sum, E_sumsquares, Mx, My, Mz, M, M_sum, M_sumsquares, M_sumfour, b, moves_accepted);

#pragma omp critical
            {
                acceptance(moves_accepted);
                Binder(M_sumsquares, M_sumfour);

                heat_cap(C, E_sum, E_sumsquares, T);
                write_C(C);
                susceptibility(chi, M_sum, M_sumsquares, T);
                write_chi(chi);

                write_M(Mx, My, Mz, M_sum);

                write_E(E_sum);

                iterate_T(T);
            }
        }
    }
    return 0;
}

//fftw instead of autocorrelation
