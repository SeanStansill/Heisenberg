//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_nearest_neighbours(int (&near_n)[L][2]);
void uniform_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]);
void random_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]);
void magnetization(double theta[L][L][L], double phi[L][L][L], double &Mx, double &My, double &Mz);
void total_energy(const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], double &E);
double local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2]);
double new_local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], const double t1[L][L][L], const double t2[L][L][L]);
void heat_cap(double &C, const double sum, const double sq_sum, const double T);
void MC_loop(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double E, double &E_sum, double &E_sumsquares);
void iterate_T(double &T);
void write_M(const double &Mx, const double &My, const double &Mz);
void write_E(const double E);
void write_C(const double C);
int random_index();
double random_ZeroTwo();
double random_decimal();
void thermalize(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T);