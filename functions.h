//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_nearest_neighbours(int (&near_n)[L][2]);
void uniform_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]);
void random_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]);
void magnetization(double theta[L][L][L], double phi[L][L][L], double &Mx, double &My, double &Mz);
void total_energy(const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], double &E);
double local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2]);
double new_local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double t1[L][L][L], double t2[L][L][L]);
void heat_cap(double &C, double E_sum, double E_sumsquares, double T);
void MC_loop(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares);
void iterate_T(double &T);
void write_M(const double Mx, const double My, const double Mz);
void write_E(double E);
void write_C(double C);
int random_index();
double random_decimal();
double random_ZeroTwo();
void thermalize(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T);
void thermalize_stat(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &sigma, double &E, double &E_sum, double &E_sumsquares);
void thermal_lattice(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T);
void MC_loopautocorr(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &E_initial, double (&zero_timeE)[steps], double (&timeE)[steps], double &ExpectationZero_T, double &Corr);
void initial_corr(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &E_initial, double (&zero_timeE)[steps], double (&timeE)[steps], double &ExpectationZero_T);
