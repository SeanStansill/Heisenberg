//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_neighbours(int near_n[N][2]);
void uniform_state(double theta[N][N][N], double phi[N][N][N]);
void random_state(double theta[N][N][N], double phi[N][N][N]);
std::vector<double> magnetization(double theta[N][N][N], double phi[N][N][N]);
double total_energy(double theta[N][N][N], double phi[N][N][N], int near_n[N][2]);
double local_energy(int i, int j, int k, double theta[N][N][N], double phi[N][N][N], int near_n[N][2]);
double new_local_energy(int i, int j, int k, double theta[N][N][N], double phi[N][N][N], int near_n[N][2], double t1[N][N][N], double t2[N][N][N]);
void write_M(std::vector<double> M);
void write_E(double E);
void write_T(double T);
void write_C(double C);
double sum_squares(double theta[N][N][N], double phi[N][N][N], int near_n[N][2]);
//double heat_cap(double E, double e_squared, double T);
double heat_cap(double theta[N][N][N], double phi[N][N][N], int near_n[N][2], double T);