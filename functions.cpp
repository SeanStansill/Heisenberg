//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include <cmath>
#include <algorithm>
#include <cstdio>

std::random_device generator;
std::mt19937 rng(generator());
std::ofstream myfile;


void write_E(const double E) {
    myfile.open("energy.txt", std::ios::app);
    myfile << E << std::endl;
    myfile.close();
}

void write_M(const double Mx, const double My, const double Mz, const double M){
    myfile.open("mag_x.txt", std::ios::app);
    myfile << Mx << std::endl;
    myfile.close();

    myfile.open("mag_y.txt", std::ios::app);
    myfile << My << std::endl;
    myfile.close();

    myfile.open("mag_z.txt", std::ios::app);
    myfile << Mz << std::endl;
    myfile.close();

    myfile.open("mag.txt", std::ios::app);
    myfile << M << std::endl;
    myfile.close();
}

void write_totalE(const double E, const int num){
    std::string name = "energy_sweep" + std::to_string(num) + ".txt";
    myfile.open(name, std::ios::app);
    myfile << E << std::endl;
    myfile.close();
}

void write_mag(const double Mx, const double My, const double Mz, const double M, const int num){
    std::string name = "mag_sweep" + std::to_string(num) + ".txt";
    myfile.open(name, std::ios::app);
    myfile << M << std::endl;
    myfile.close();
}

void write_N(const int N){
    myfile.open("N.txt", std::ios::app);
    myfile << N << std::endl;
    myfile.close();
}

void write_Eloop(const double E, const std::string name){

    myfile.open(name, std::ios::app);
    myfile << E << std::endl;
    myfile.close();
}

void get_nearest_neighbours(int (&near_n)[L][2]){
    for (int l = 0; l < L; l++) {
        near_n[l][0] = ((l + 1) % L);
        near_n[l][1] = ((l - 1 + L) % L);
    }
}

int random_index(){
    std::uniform_int_distribution<int> uniform_distribution;
    int index;
    index = uniform_distribution(rng);
    return index;
}

double random_decimal(){
    std::uniform_real_distribution<double> decimal(-1.0*pi, 1.0*pi);
    double dec = decimal(rng);
    return dec;
}

double random_ZeroTwo(){
    std::uniform_real_distribution<double> ZeroTwo(-2.0*pi, 2.0*pi);
    double angle = ZeroTwo(rng);
    return angle;

}

double random_uniform(){
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    double num = uniform_distribution(rng);
    return num;

}

void random_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]){
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                double ang1 =  random_ZeroTwo(); //random_decimal() * 2.0 * pi;
                double ang2 =  random_decimal(); //acos((2.0 * random_decimal())-1);
                theta[l][m][n] = ang1; // This corrects for overcounting at the poles of the sphere. See: http://mathworld.wolfram.com/SpherePointPicking.html
                phi[l][m][n] = ang2;
            }
        }
    }
}

void uniform_spin_state(double (&theta)[L][L][L], double (&phi)[L][L][L]){
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                theta[l][m][n] = 0.0;
                phi[l][m][n] = pi/2;
            }
        }
    }
}

void magnetization(const double theta[L][L][L], const double phi[L][L][L], double &Mx, double &My, double &Mz, double &M){
    Mx = 0.0;
    My = 0.0;
    Mz = 0.0;
    M = 0.0;
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                Mx += cos(theta[l][m][n]) * sin(phi[l][m][n]);
                My += sin(theta[l][m][n]) * sin(phi[l][m][n]);
                Mz += cos(phi[l][m][n]);
            }
        }
    }
    Mx = (Mx/N);
    My = (My/N);
    Mz = (Mz/N);
    M = sqrt((Mx*Mx) + (My*My) + (Mz*Mz));
}

void local_M(const double theta, const double phi, const double theta_trial, const double phi_trial, double &dSx, double &dSy, double &dSz, double &dS){
    dSx = 0.0;
    dSy = 0.0;
    dSz = 0.0;
    dS = 0.0;
    dSx = (cos(theta_trial) * sin(phi_trial)) - (cos(theta) * sin(phi));
    dSy = (sin(theta_trial) * sin(phi_trial)) - (sin(theta) * sin(phi));
    dSz = cos(phi_trial) - cos(phi);
    dS = dSx + dSy + dSz;
}

void susceptibility(double &chi, const double M_sum, const double M_sumsquares, const double T){
    chi = (M_sumsquares-(M_sum*M_sum))/(T*T);
}


double local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2]){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]);
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]);
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]);
        e += sin(theta[i][j][k]) * sin(phi[i][j][k]) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += sin(theta[i][j][k]) * sin(phi[i][j][k]) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += sin(theta[i][j][k]) * sin(phi[i][j][k]) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += cos(phi[i][j][k]) * cos(phi[near_n[i][o]][j][k]);
        e += cos(phi[i][j][k]) * cos(phi[i][near_n[j][o]][k]);
        e += cos(phi[i][j][k]) * cos(phi[i][j][near_n[k][o]]);
    }
    return e * -J;
}
double new_local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], const double theta_trial[L][L][L] = {}, const double phi_trial[L][L][L] = {}){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += cos(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]);
        e += cos(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]);
        e += cos(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]);
        e += sin(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += sin(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += sin(theta_trial[i][j][k]) * sin(phi_trial[i][j][k]) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += cos(phi_trial[i][j][k]) * cos(phi[near_n[i][o]][j][k]);
        e += cos(phi_trial[i][j][k]) * cos(phi[i][near_n[j][o]][k]);
        e += cos(phi_trial[i][j][k]) * cos(phi[i][j][near_n[k][o]]);
    }
    return e * -J;
}

void total_energy(const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], double &E) {
    E = 0.0;
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {

                    E += local_energy(l, m, n, theta, phi, near_n);
            }
        }
    }
    //E = E/N;
}

void heat_cap(double &C, const double E_sum, const double E_sumsquares, const double T){
    C = (E_sumsquares - (E_sum*E_sum))/(T*T);
}

void Binder(double &M_sumsquares, double &M_sumfour){

    double U = 1.0 - (M_sumfour / (3 * M_sumsquares * M_sumsquares));

    myfile.open("binder.txt", std::ios::app);
    myfile << U << std::endl;
    myfile.close();
}


void MC_loop(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &Mx, double &My, double &Mz, double &M, double &M_sum, double &M_sumsquares, double &M_sumfour, const int b, double &Sx, double &Sy, double &Sz, double &S) {
    E_sum = 0.0;
    E_sumsquares = 0.0;
    M_sum = 0.0;
    M_sumsquares = 0.0;
    M_sumfour = 0.0;
    int i, j, k, counter = 0;
    double E1, E2, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {}, moves_accepted = 0.0;
    double dE = 0.0;
    for (int a = 0; a < nsteps; a++) {
        i = random_index() % L;
        j = random_index() % L;
        k = random_index() % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = 2.0*pi*random_uniform(); //theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
        phi_trial[i][j][k] = std::acos((2.0*random_uniform()) - 1); //phi[i][j][k] + (0.05 * random_decimal() * pi);


        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if (random_uniform() <= exp(std::min(0.0, -dE / T))) {
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
            moves_accepted += 1.0;
            //local_M(theta[i][j][k], phi[i][j][k], Sx, Sy, Sz, S);
            //E = ((E*N) + dE)/N;
            //Mx = ((Mx*N) + Sx)/N;
            //My = ((My*N) + Sy)/N;
            //Mz = ((Mz*N) + Sz)/N;
            //M = sqrt((Mx*Mx) + (My*My) + (Mz*Mz));
        }


        if (a % 1 == 0) {
            total_energy(theta, phi, near_n, E);  // Will become obsolete once local updates are implemented
            E_sum += E;
            E_sumsquares += (E * E);
            magnetization(theta, phi, Mx, My, Mz, M); // Will become obsolete once local updates are implemented
            M_sum += M;
            M_sumsquares += (M * M);
            M_sumfour += (M * M * M * M);
            counter++;
        }
    }
    E_sum = E_sum / counter;
    E_sumsquares = E_sumsquares / counter;
    M_sum = M_sum / counter;
    M_sumsquares = M_sumsquares / counter;
    M_sumfour = M_sumfour / counter;
    myfile.open("acceptance.txt", std::ios::app);
    myfile << moves_accepted/nsteps << std::endl;
    myfile.close();
    Binder(M_sumsquares, M_sumfour);
//}
}

void MC_typewriter(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &Mx, double &My, double &Mz, double &M, double &M_sum, double &M_sumsquares, double &M_sumfour, const int b, double &dSx, double &dSy, double &dSz, double &dS) {
    E_sum = 0.0;
    E_sumsquares = 0.0;
    M_sum = 0.0;
    M_sumsquares = 0.0;
    M_sumfour = 0.0;
    int i, j, k, counter = 0;
    double E1, E2, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {}, moves_accepted = 0.0;
    double dE = 0.0;
    for (int a = 0; a < sweeps; a++) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < L; k++) {

                    E1 = local_energy(i, j, k, theta, phi, near_n);

                    theta_trial[i][j][k] =
                            2.0 * pi * random_uniform(); //theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
                    phi_trial[i][j][k] = std::acos(
                            (2.0 * random_uniform()) - 1); //phi[i][j][k] + (0.05 * random_decimal() * pi);


                    E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

                    dE = E2 - E1;
                    if (random_uniform() <= exp(std::min(0.0, -dE / T))) {
                        //local_M(theta[i][j][k], phi[i][j][k], theta_trial[i][j][k], phi_trial[i][j][k], dSx, dSy, dSz, dS);
                        //E = E + dE;
                        //Mx = ((Mx*N) + dSx)/N;
                        //My = ((My*N) + dSy)/N;
                        //Mz = ((Mz*N) + dSz)/N;
                        //M = sqrt((Mx*Mx) + (My*My) + (Mz*Mz));
                        theta[i][j][k] = theta_trial[i][j][k];
                        phi[i][j][k] = phi_trial[i][j][k];
                        moves_accepted += 1.0;
                    }


                    if (a % 5 == 0) {
                        total_energy(theta, phi, near_n, E);  // Will become obsolete once local updates are implemented
                        E_sum += E;
                        E_sumsquares += (E * E);
                        magnetization(theta, phi, Mx, My, Mz,M); // Will become obsolete once local updates are implemented
                        M_sum += M;
                        M_sumsquares += (M * M);
                        M_sumfour += (M * M * M * M);
                        counter++;
                    }
                }
            }
        }
    }
    E_sum = E_sum / counter;
    E_sumsquares = E_sumsquares / counter;
    M_sum = M_sum / counter;
    M_sumsquares = M_sumsquares / counter;
    M_sumfour = M_sumfour / counter;
    myfile.open("acceptance.txt", std::ios::app);
    myfile << moves_accepted/nsteps << std::endl;
    myfile.close();
    Binder(M_sumsquares, M_sumfour);
//}
}

void thermalize_typewriter(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T){
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < thermalsteps; a++) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < L; k++) {

                    E1 = local_energy(i, j, k, theta, phi, near_n);

                    theta_trial[i][j][k] = 2.0 * pi * random_uniform();
                    //theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
                    phi_trial[i][j][k] = std::acos((2.0 * random_uniform()) - 1);
                    //phi[i][j][k] + (0.05 * random_decimal() * pi);



                    E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

                    dE = E2 - E1;
                    if (dE < 0) {
                        theta[i][j][k] = theta_trial[i][j][k];
                        phi[i][j][k] = phi_trial[i][j][k];
                    } else if (random_uniform() < exp(std::min(0.0, -dE / T))) {
                        theta[i][j][k] = theta_trial[i][j][k];
                        phi[i][j][k] = phi_trial[i][j][k];
                    }
                }
            }
        }
    }
}

void iterate_T(double &T){
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
    T += T_step;
}

void write_C(const double C){
    myfile.open("heat_cap.txt", std::ios::app);
    myfile << C << std::endl;
    myfile.close();
}

void write_chi(const double chi){
    myfile.open("susceptibility.txt", std::ios::app);
    myfile << chi << std::endl;
    myfile.close();
}

void write_Corr(const double &Corr){
    myfile.open("correlation.txt", std::ios::app);
    myfile << Corr << std::endl;
    myfile.close();
}