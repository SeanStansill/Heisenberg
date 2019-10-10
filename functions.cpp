//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include <cmath>
#include <algorithm>

std::random_device generator;
std::mt19937 rng(generator());
std::ofstream myfile;


void write_E(const double E) {
    myfile.open("energy.txt", std::ios::app);
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
    std::uniform_real_distribution<double> decimal(0.0, 1.0*pi);
    double dec = decimal(rng);
    return dec;
}

double random_ZeroTwo(){
    std::uniform_real_distribution<double> ZeroTwo(0.0, 2*pi);
    double angle = ZeroTwo(rng);
    return angle;

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

void magnetization(double theta[L][L][L], double phi[L][L][L], double &Mx, double &My, double &Mz, double &M){
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

void susceptibility(double &chi, const double M_sum, const double M_sumsquares, const double T){
    chi = (M_sumsquares-M_sum)/(T*T);
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
    E = E/N;
}

void heat_cap(double &C, const double E_sum, const double E_sumsquares, const double T){
    C = (E_sumsquares - (E_sum*E_sum))/(T*T);
}

void thermalize(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T){
    int i, j, k;
    double Et;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < thermalsteps; a++) {
        i = random_index() % L;
        j = random_index() % L;
        k = random_index() % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
        phi_trial[i][j][k] = phi[i][j][k] + (0.05 * random_decimal() * pi);

        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if(random_decimal() <= exp(std::min(0.0, -dE/T))){
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }
        /*total_energy(theta, phi, near_n, Et);
        myfile.open("trial.txt", std::ios::app);
        myfile << Et << std::endl;
        myfile.close();*/
    }
}

void MC_loop(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &Mx, double &My, double &Mz, double &M, double &M_sum, double &M_sumsquares){
    E_sum = 0.0;
    E_sumsquares = 0.0;
    M_sum = 0.0;
    M_sumsquares = 0.0;
    int i, j, k, counter = 0;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < nsteps; a++) {
        i = random_index() % L;
        j = random_index() % L;
        k = random_index() % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
        phi_trial[i][j][k] = phi[i][j][k] + (0.05 * random_decimal() * pi);



        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if(random_decimal() <= exp(std::min(0.0, -dE/T))){
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }

        if(a % 5 == 0) {
            total_energy(theta, phi, near_n, E);
            E_sum += E;
            E_sumsquares += (E * E);
            magnetization(theta, phi, Mx, My, Mz, M);
            M_sum += M;
            M_sumsquares += (M*M);
            counter++;
        }
    }
    E_sum = E_sum/counter;
    E_sumsquares = E_sumsquares/counter;
    M_sum = M_sum/counter;
    M_sumsquares = M_sum/counter;
}

void iterate_T(double &T){
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
    T += T_step;
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

void thermal_lattice(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T) {
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < latticesweeps; a++) {
        for (int l = 0; l < L; l++) {
            for (int m = 0; m < L; m++) {
                for (int n = 0; n < L; n++) {

                    E1 = local_energy(l, m, n, theta, phi, near_n);

                    theta_trial[l][m][n] = theta[l][m][n] + (0.05 * random_ZeroTwo() * pi);
                    phi_trial[l][m][n] = phi[l][m][n] + (0.05 * random_decimal() * pi);

                    E2 = new_local_energy(l, m, n, theta, phi, near_n, theta_trial, phi_trial);

                    dE = E2 - E1;
                    if (random_decimal() <= exp(std::min(0.0, (((-J) / (k_B * T)) * dE)))) {
                        theta[l][m][n] = theta_trial[l][m][n];
                        phi[l][m][n] = phi_trial[l][m][n];
                    }
                }
            }
        }
    }
}

void MC_loopautocorr(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &E_initial, double (&zero_timeE)[steps], double (&timeE)[steps], double &ExpectationZero_T, double &Corr) {
    E_sum = 0.0;
    E_sumsquares = 0.0;
    Corr = 10.0;
    total_energy(theta, phi, near_n, E);
    int i, j, k;

    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    while (Corr > (1.0 /(2*exp1))) {
        for (int a = 0; a < steps; a++) {
            i = random_index() % L;
            j = random_index() % L;
            k = random_index() % L;

            E1 = local_energy(i, j, k, theta, phi, near_n);

            theta_trial[i][j][k] = theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
            phi_trial[i][j][k] = phi[i][j][k] + (0.05 * random_decimal() * pi);


            E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

            dE = E2 - E1;
            if (random_decimal() <= exp(std::min(0.0, -dE/T))){
                theta[i][j][k] = theta_trial[i][j][k];
                phi[i][j][k] = phi_trial[i][j][k];
            }

            total_energy(theta, phi, near_n, E);
            E_sum += E;
            E_sumsquares += (E*E);
            timeE[a] = E;
        }
        ExpectationZero_T = 0.0;
        for(int l = 0; l < steps; l++){
        ExpectationZero_T += (timeE[l]*zero_timeE[l]);
        }
        Corr = ((ExpectationZero_T) - (E_initial*E_sum))/(E_sumsquares - (E_sum*E_sum));
        write_Corr(Corr);
    }
}

void initial_corr(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares, double &E_initial, double (&zero_timeE)[steps], double (&timeE)[steps], double &ExpectationZero_T) {
    E_sum = 0.0;
    E_sumsquares = 0.0;
    int i, j, k;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
        for (int a = 0; a < steps; a++) {
            i = random_index() % L;
            j = random_index() % L;
            k = random_index() % L;

            E1 = local_energy(i, j, k, theta, phi, near_n);

            theta_trial[i][j][k] = theta[i][j][k] + (0.05 * random_ZeroTwo() * pi);
            phi_trial[i][j][k] = phi[i][j][k] + (0.05 * random_decimal() * pi);


            E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

            dE = E2 - E1;
            if (random_decimal() <= exp(std::min(0.0, -dE/T))){
                theta[i][j][k] = theta_trial[i][j][k];
                phi[i][j][k] = phi_trial[i][j][k];
            }

            total_energy(theta, phi, near_n, E);
            E_sum += E;
            E_sumsquares += (E*E);
            timeE[a] = E;
        }
        ExpectationZero_T = 0.0;
        for(int l = 0; l < steps; l++){
            ExpectationZero_T += (timeE[l]*zero_timeE[l]);
        }
        E_initial = E_sum;
        for(int l = 0; l < steps; l++) {
            zero_timeE[l] = timeE[l];
        }
    }
