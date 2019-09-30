//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include <fstream>


void get_neighbours(int near_n[N][2]){
    for (int l = 0; l < N; l++) {
        near_n[l][0] = ((l + 1) % N);
        near_n[l][1] = ((l - 1 + N) % N);
    }
}


void random_state(double theta[N][N][N], double phi[N][N][N]){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> ang(0.0, 2.0);
    for (int l = 0; l < N; l++) {
        for (int m = 0; m < N; m++) {
            for (int n = 0; n < N; n++) {
                double ang1 = ang(mt);
                double ang2 = ang(mt);
                theta[l][m][n] = ang1*pi;
                phi[l][m][n] = ang2*pi;
            }
        }
    }
}

void uniform_state(double theta[N][N][N], double phi[N][N][N]){
    for (int l = 0; l < N; l++) {
        for (int m = 0; m < N; m++) {
            for (int n = 0; n < N; n++) {
                theta[l][m][n] = 0.0;
                phi[l][m][n] = pi/2;
            }
        }
    }
}

std::vector<double> magnetization(double theta[N][N][N], double phi[N][N][N]) {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    std::vector<double> sum;

    for (int l = 0; l < N; ++l) {
        for (int m = 0; m < N; ++m) {
            for (int n = 0; n < N; ++n) {
                x += cos(theta[l][m][n]) * sin(phi[l][m][n]);
                y += sin(theta[l][m][n]) * sin(phi[l][m][n]);
                z += cos(phi[l][m][n]);

            }
        }
    }
    x *= (1/L);
    y *= (1/L);
    z *= (1/L);
    sum.push_back(x);
    sum.push_back(y);
    sum.push_back(z);
    return sum;
}

double total_energy(double theta[N][N][N], double phi[N][N][N], int near_n[N][2]){
    double E_t = 0.0;
    for (int l = 0; l < N; l++) {
        for (int m = 0; m < N; m++) {
            for (int n = 0; n < N; n++) {
                for(int o = 0; o < 2; o++) {

                    E_t += (cos(theta[l][m][n]) * sin(phi[l][m][n])) * (cos(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]));
                    E_t += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]);
                    E_t += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]);
                    E_t += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]));
                    E_t += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]));
                    E_t += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]));
                    E_t += cos(phi[l][m][n]) * cos(phi[near_n[l][o]][m][n]);
                    E_t += cos(phi[l][m][n]) * cos(phi[l][near_n[m][o]][n]);
                    E_t += cos(phi[l][m][n]) * cos(phi[l][m][near_n[n][o]]);
                }
            }
        }
    }
    E_t *= (1/L)*(-J);
    return E_t;
}

double local_energy(int i, int j, int k, double theta[N][N][N], double phi[N][N][N], int near_n[N][2]){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]);
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]);
        e += cos(theta[i][j][k]) * sin(phi[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]);
        e += (sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += (sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += (sin(theta[i][j][k]) * sin(phi[i][j][k])) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += cos(phi[i][j][k]) * cos(phi[near_n[i][o]][j][k]);
        e += cos(phi[i][j][k]) * cos(phi[i][near_n[j][o]][k]);
        e += cos(phi[i][j][k]) * cos(phi[i][j][near_n[k][o]]);
    }
    e *= (1 / L) * (-J);
    return e;
}
double new_local_energy(int i, int j, int k, double theta[N][N][N], double phi[N][N][N], int near_n[N][2], double t1[N][N][N] = {}, double t2[N][N][N] = {}){
    double e = 0.0;
    for(int o = 0; o<2; o++) {
        e += cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]);
        e += cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]);
        e += cos(t1[i][j][k]) * sin(t2[i][j][k]) * cos(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]);
        e += (sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[near_n[i][o]][j][k]) * sin(phi[near_n[i][o]][j][k]));
        e += (sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[i][near_n[j][o]][k]) * sin(phi[i][near_n[j][o]][k]));
        e += (sin(t1[i][j][k]) * sin(t2[i][j][k])) * (sin(theta[i][j][near_n[k][o]]) * sin(phi[i][j][near_n[k][o]]));
        e += cos(t2[i][j][k]) * cos(phi[near_n[i][o]][j][k]);
        e += cos(t2[i][j][k]) * cos(phi[i][near_n[j][o]][k]);
        e += cos(t2[i][j][k]) * cos(phi[i][j][near_n[k][o]]);
    }
    e *= (1 / L) * (-J);
    return e;
}

double sum_squares(double theta[N][N][N], double phi[N][N][N], int near_n[N][2]){
    double e_squared = 0.0;
    for (int l = 0; l < N; l++) {
        for (int m = 0; m < N; m++) {
            for (int n = 0; n < N; n++) {
                for(int o = 0; o < 2; o++) {
                    double E_1 = 0.0;
                    double E_2 = 0.0;
                    double E_3 = 0.0;

                    E_1 += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]);
                    E_1 += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[near_n[l][o]][m][n]) * sin(phi[near_n[l][o]][m][n]));
                    E_1 += cos(phi[l][m][n]) * cos(phi[near_n[l][o]][m][n]);
                    E_1 *= (1/L)*(-J);

                    E_2 += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]);
                    E_2 += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][near_n[m][o]][n]) * sin(phi[l][near_n[m][o]][n]));
                    E_2 += cos(phi[l][m][n]) * cos(phi[l][near_n[m][o]][n]);
                    E_2 *= (1/L)*(-J);

                    E_3 += cos(theta[l][m][n]) * sin(phi[l][m][n]) * cos(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]);
                    E_3 += (sin(theta[l][m][n]) * sin(phi[l][m][n])) * (sin(theta[l][m][near_n[n][o]]) * sin(phi[l][m][near_n[n][o]]));
                    E_3 += cos(phi[l][m][n]) * cos(phi[l][m][near_n[n][o]]);
                    E_3 *= (1/L)*(-J);

                    e_squared += (E_1*E_1) + (E_2*E_2) + (E_3*E_3);
                }
            }
        }
    }
    return e_squared;
}

double heat_cap(double theta[N][N][N], double phi[N][N][N], int near_n[N][2], double T){  //(double E, double e_squared, double T){
    double E = total_energy(theta, phi, near_n);
    double e_squared = sum_squares(theta, phi, near_n);
    double C = (e_squared - (E*E))/(k_B*T*T);
    return C;
}

void write_M(std::vector<double> M){
    std::ofstream myfile;

    myfile.open("mag_x.txt", std::ios::app);
    myfile << M[0];
    myfile << std::endl;
    myfile.close();

    myfile.open("mag_y.txt", std::ios::app);
    myfile << M[1];
    myfile << std::endl;
    myfile.close();

    myfile.open("mag_z.txt", std::ios::app);
    myfile << M[2];
    myfile << std::endl;
    myfile.close();
}

void write_E(double E){
    std::ofstream myfile;
    myfile.open("energy.txt", std::ios::app);
    myfile << E << std::endl;
    myfile.close();
}

void write_C(double C){
    std::ofstream myfile;
    myfile.open("heat_cap.txt", std::ios::app);
    myfile << C << std::endl;
    myfile.close();
}

void write_T(double T){
    std::ofstream myfile;
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
}
