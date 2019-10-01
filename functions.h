//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_neighbours(int near_n[L][2]);
void uniform_state(double theta[L][L][L], double phi[L][L][L]);
void random_state(double theta[L][L][L], double phi[L][L][L]);
std::vector<double> magnetization(double theta[L][L][L], double phi[L][L][L]);
double total_energy(double theta[L][L][L], double phi[L][L][L], int near_n[L][2]);
double local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2]);
double new_local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double t1[L][L][L], double t2[L][L][L]);
void write_M(std::vector<double> M);
void write_E(double E);
void write_T(double T);
void write_C(double C);
double heat_cap(double sum, double sq_sum, double T);
//double heat_cap(double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double T);