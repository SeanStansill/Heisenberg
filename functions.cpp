//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include <fstream>
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

void magnetization(double theta[L][L][L], double phi[L][L][L], double &Mx, double &My, double &Mz){
    Mx = 0.0;
    My = 0.0;
    Mz = 0.0;
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
}

void total_energy(const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], double &E) {
    E = 0.0;
    for (int l = 0; l < L; l++) {
        for (int m = 0; m < L; m++) {
            for (int n = 0; n < L; n++) {
                for (int o = 0; o < 2; o++) {

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
    E *= ((-J) / N);
}

double local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2]){
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
double new_local_energy(const int i, const int j, const int k, const double theta[L][L][L], const double phi[L][L][L], const int near_n[L][2], const double t1[L][L][L] = {}, const double t2[L][L][L] = {}){
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

void heat_cap(double &C, const double E_sum, const double E_sumsquares, const double T){
    C = (E_sumsquares - (E_sum*E_sum))/(k_B*T*T);
}

void thermalize(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T){
    int i, j, k;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < thermalsteps; a++) {
        i = random_index() % L;
        j = random_index() % L;
        k = random_index() % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = random_ZeroTwo() * pi;
        phi_trial[i][j][k] = random_ZeroTwo() * pi;

        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if(random_decimal() <= exp(std::min(0.0, (((-J)/(k_B*T)) * dE)))){
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }
    }
}

void MC_loop(double (&theta)[L][L][L], double (&phi)[L][L][L], const int near_n[L][2], const double T, double &E, double &E_sum, double &E_sumsquares){
    E_sum = 0.0;
    E_sumsquares = 0.0;
    int i, j, k, counter = 0;
    double E1, E2, dE, theta_trial[L][L][L] = {}, phi_trial[L][L][L] = {};
    for (int a = 0; a < nsteps; a++) {
        i = random_index() % L;
        j = random_index() % L;
        k = random_index() % L;

        E1 = local_energy(i, j, k, theta, phi, near_n);

        theta_trial[i][j][k] = theta[i][j][k] + (0.01 * random_ZeroTwo() * pi);
        phi_trial[i][j][k] = phi[i][j][k] + (0.01 * random_ZeroTwo() * pi);

        //theta_trial[i][j][k] = random_ZeroTwo() * pi;
        //phi_trial[i][j][k] = random_ZeroTwo() * pi;

        E2 = new_local_energy(i, j, k, theta, phi, near_n, theta_trial, phi_trial);

        dE = E2 - E1;
        if(random_decimal() <= exp(std::min(0.0, (((-J)/(k_B*T)) * dE)))){
            theta[i][j][k] = theta_trial[i][j][k];
            phi[i][j][k] = phi_trial[i][j][k];
        }

        //total_energy(theta, phi, near_n, E);
        //write_E(E);

        if(a % 10 == 0) {
            total_energy(theta, phi, near_n, E);
            E_sum += E;
            E_sumsquares += (E * E);
            counter++;
        }
    }
    //E_sum = E_sum/nsteps;
    //E_sumsquares = E_sumsquares/nsteps;
    E_sum = E_sum/counter;
    E_sumsquares = E_sumsquares/counter;
}

void iterate_T(double &T){
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
    T += T_step;
}

void write_M(const double Mx, const double My, const double Mz){
    myfile.open("mag_x.txt", std::ios::app);
    myfile << Mx << std::endl;
    myfile.close();

    myfile.open("mag_y.txt", std::ios::app);
    myfile << My << std::endl;
    myfile.close();

    myfile.open("mag_z.txt", std::ios::app);
    myfile << Mz << std::endl;
    myfile.close();
}




void write_C(const double C){
    myfile.open("heat_cap.txt", std::ios::app);
    myfile << C << std::endl;
    myfile.close();
}
