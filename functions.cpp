//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include <fstream>
#include <cmath>
#include <algorithm>
#include <random>


void get_neighbours(int near_n[L][2]){
    for (int l = 0; l < L; l++) {
        near_n[l][0] = ((l + 1) % L);
        near_n[l][1] = ((l - 1 + L) % L);
    }
}


void random_state(double theta[L][L][L], double phi[L][L][L]){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> ang(0.0, 2.0);
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                double ang1 = ang(mt);
                double ang2 = ang(mt);
                theta[l][m][n] = ang1*pi; // This corrects for overcounting at the poles of the sphere. See: http://mathworld.wolfram.com/SpherePointPicking.html
                phi[l][m][n] = ang2*pi;
            }
        }
    }
}

void uniform_state(double theta[L][L][L], double phi[L][L][L]){
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                theta[l][m][n] = 0.0;
                phi[l][m][n] = pi/2;
            }
        }
    }
}

void magnetization(double theta[L][L][L], double phi[L][L][L]) {
    double Mx = 0.0;
    double My = 0.0;
    double Mz = 0.0;

    for (int l = 0; l < L; ++l) {
        for (int m = 0; m < L; ++m) {
            for (int n = 0; n < L; ++n) {
                Mx += cos(theta[l][m][n]) * sin(phi[l][m][n]);
                My += sin(theta[l][m][n]) * sin(phi[l][m][n]);
                Mz += cos(phi[l][m][n]);

            }
        }
    }
    Mx *= (1/N);
    My *= (1/N);
    Mz *= (1/N);

    std::ofstream myfile;

    myfile.open("mag_x.txt", std::ios::app);
    myfile << Mx;
    myfile << std::endl;
    myfile.close();

    myfile.open("mag_y.txt", std::ios::app);
    myfile << My;
    myfile << std::endl;
    myfile.close();

    myfile.open("mag_z.txt", std::ios::app);
    myfile << Mz;
    myfile << std::endl;
    myfile.close();
}

void total_energy(double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double E){
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                for(int o = 0; o < 2; o++) {

                    E += (cos(theta[l][m][n]) * sin(phi[l][m][n])) * (cos(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]));
                    E += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]);
                    E += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]);
                    E += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]));
                    E += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]));
                    E += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]));
                    E += cos(phi[l][m][n]) * cos(phi[near_n[l][o]][m][n]);
                    E += cos(phi[l][m][n]) * cos(phi[l][near_n[m][o]][n]);
                    E += cos(phi[l][m][n]) * cos(phi[l][m][near_n[n][o]]);
                }
            }
        }
    }
    E *= ((-J)/N);
    std::ofstream myfile;
    myfile.open("energy.txt", std::ios::app);
    myfile << E << std::endl;
    myfile.close();
}

double local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2]){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += ((-J)/L)*(cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += ((-J)/L)*(cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += ((-J)/L)*(cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += ((-J)/L)*((sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k])));
        e += ((-J)/L)*((sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k])));
        e += ((-J)/L)*((sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]])));
        e += ((-J)/L)*(cos(phi[i][j][k]) * cos(phi[near_n[i][o]][j][k]));
        e += ((-J)/L)*(cos(phi[i][j][k]) * cos(phi[i][near_n[j][o]][k]));
        e += ((-J)/L)*(cos(phi[i][j][k]) * cos(phi[i][j][near_n[k][o]]));
    }
    return e;
}
double new_local_energy(int i, int j, int k, double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double t1[L][L][L] = {}, double t2[L][L][L] = {}){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += ((-J)/L)*(cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += ((-J)/L)*(cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += ((-J)/L)*(cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += ((-J)/L)*((sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k])));
        e += ((-J)/L)*((sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k])));
        e += ((-J)/L)*((sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]])));
        e += ((-J)/L)*(cos(t2[i][j][k]) * cos(phi[near_n[i][o]][j][k]));
        e += ((-J)/L)*(cos(t2[i][j][k]) * cos(phi[i][near_n[j][o]][k]));
        e += ((-J)/L)*(cos(t2[i][j][k]) * cos(phi[i][j][near_n[k][o]]));
    }
    return e;
}

void heat_cap(double E_sum, double E_sumsquares, double T){
    double C = (E_sumsquares - (E_sum*E_sum))/(k_B*T*T);
    std::ofstream myfile;
    myfile.open("heat_cap.txt", std::ios::app);
    myfile << C << std::endl;
    myfile.close();
}

void MC_steps(double theta[L][L][L], double phi[L][L][L], int near_n[L][2], double T, double E, double E_sum, double E_sumsquares){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> uniform_distribution;
    std::uniform_real_distribution<double> decimal(0.0, 1.0);
    int i, j, k;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < nsteps; a++) {
        i = uniform_distribution(mt) % L;
        j = uniform_distribution(mt) % L;
        k = uniform_distribution(mt) % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = 2*pi*decimal(mt);
        phi_trial[i][j][k] = acos((2*decimal(mt))-1);

        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if(decimal(mt) <= exp(std::min(0.0, (((-J)/(k_B*T)) * dE)))){
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }

        total_energy(theta, phi, near_n, E);
        E_sum += E;
        E_sumsquares += (E*E);
    }
    sum = sum/nsteps;
    E_sum = E_sumsquares/nsteps;
}

void iterate_T(double T){
    std::ofstream myfile;
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
    T += T_step;
}