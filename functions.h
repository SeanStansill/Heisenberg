//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_nearest_neighbours(int (&up)[N], int (&down)[N], int (&left)[N], int (&right)[N], int (&forwards)[N], int (&backwards)[N]);
void uniform_spin_state(double (&theta)[N], double (&phi)[N]);
void random_spin_state(double (&theta)[N], double (&phi)[N]);
void magnetization(const double theta[N], const double phi[N], double &Mx, double &My, double &Mz, double &M);
void single_spin_M(const double &theta, const double &phi, double &Mx, double &My, double &Mz, double);
double local_energy(const int i, const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N]);
double single_spin_energy(const int i, const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N]);
void total_energy(const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], double &E);
void heat_cap(double &C, double E_sum, double E_sumsquares, const double T);
void Binder(const double &M_sumsquares, const double &M_sumfour);
void MC_parallel(double (&theta)[N], double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], const double T, double &E, double &E_sum, double &E_sumsquares, double &Mx, double &My, double &Mz, double &M, double &M_sum, double &M_sumsquares, double &M_sumfour, const int b, double &moves_accepted);
void iterate_T(const double &T);
void write_M(const double Mx, const double My, const double Mz, const double M);
void write_E(const double E);
void write_C(const double C);
int random_index();
double random_decimal();
double random_ZeroTwo();
void susceptibility(double &chi, const double M_sum, const double M_sumsquares, const double T);
void write_chi(const double chi);
void write_N(const int N);
static inline double random_uniform();
void thermalize_parallel(double (&theta)[N], double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], const double T);
void acceptance(const double moves_accepted);
