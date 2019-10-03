
#include "constants.h"
#include "functions.h"

int main() {

    double theta[L][L][L], phi[L][L][L], E, T = T_init, E_sum, E_sumsquares, C, Mx, My, Mz;
    int near_n[L][2];

    get_nearest_neighbours(near_n);

    //random_spin_state(theta, phi);
    uniform_spin_state(theta, phi);

    magnetization(theta, phi, Mx, My, Mz);
    write_M(Mx, My, Mz);

    total_energy(theta, phi, near_n, E);
    write_E(E);


    for(int b = 0; b < T_intervals; b++) {
        thermalize(theta, phi, near_n, T);
        MC_loop(theta, phi, near_n, T, E, E_sum, E_sumsquares);

        heat_cap(C, E_sum, E_sumsquares, T);
        write_C(C);

        magnetization(theta, phi, Mx, My, Mz);

        total_energy(theta, phi, near_n, E);
        write_E(E);

        iterate_T(T);
    }
}
