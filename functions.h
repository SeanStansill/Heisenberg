//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

void get_neighbours(int near_n[L][2]);
void uniform_state(double theta[L][L][L], double phi[L][L][L]);
void random_state(double theta[L][L][L], double phi[L][L][L]);
void magnetization(double theta[L][L][L], double phi[L][L][L]);
void total_energy(double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double E);
double local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2]);
double new_local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double t1[L][L][L], double t2[L][L][L]);
void heat_cap(double sum, double sq_sum, double T);
void MC_steps(double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double T, double E, double E_sum, double E_sumsquares);
void iterate_T(double T);
