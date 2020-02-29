//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <random>
#include"constants.h"
#include "functions.h"

#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>

std::ofstream myfile;

void acceptance(const double moves_accepted){
    myfile.open("acceptance.txt", std::ios::app);
    myfile << moves_accepted/(sweeps*N) << std::endl;
    myfile.close();
}

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

void write_N(const int N){
    myfile.open("N.txt", std::ios::app);
    myfile << N << std::endl;
    myfile.close();
}

void get_nearest_neighbours(int (&up)[N], int (&down)[N], int (&left)[N], int (&right)[N], int (&forwards)[N], int (&backwards)[N]){
    myfile.open("neighbours.txt", std::ios::app);
    myfile << "Index" << "    " << "Up" << "    " << "Down" << "    " << "Left" << "    " << "Right" << "    " << "Backwards" << "    " << "Forwards" << std::endl;
    myfile.close();
    for (int i = 0; i < N; i++) {
        if(i % (L*L) > L*(L-1) - 1){up[i] = i-(L*(L-1));}
        else{up[i] = i + L;}

        if((i % (L*L)) < L){down[i] = i + (L*(L-1));}
        else{down[i] = i - L;}

        if(i % L == 0){left[i] = i + (L-1);}
        else{left[i] = i - 1;}

        if((i + 1) % L == 0){right[i] = i - (L-1);}
        else{right[i] = i + 1;}

        forwards[i] = ((i - (L*L)) + N) % N;

        backwards[i] = (i + (L*L)) % N;

        myfile.open("neighbours.txt", std::ios::app);
        myfile << i << "    " << up[i] << "    " << down[i] << "    " << left[i] << "    " << right[i] << "    " << backwards[i] << "    " << forwards[i] << std::endl;
        myfile.close();
    }
}

int random_index(){
    static thread_local std::random_device generator;
    static thread_local std::mt19937 rng(generator());
    std::uniform_int_distribution<int> uniform_distribution;
    int index;
    index = uniform_distribution(rng);
    return index;
}

double random_decimal(){
    static thread_local std::random_device generator;
    static thread_local std::mt19937 rng(generator());
    std::uniform_real_distribution<double> decimal(-1.0*pi, 1.0*pi);
    double dec = decimal(rng);
    return dec;
}

double random_ZeroTwo(){
    static thread_local std::random_device generator;
    static thread_local std::mt19937 rng(generator());
    std::uniform_real_distribution<double> ZeroTwo(-2.0*pi, 2.0*pi);
    double angle = ZeroTwo(rng);
    return angle;

}

double random_uniform(){
    static thread_local std::random_device generator;
    static thread_local std::mt19937 rng(generator());
    std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);
    double num = uniform_distribution(rng);
    return num;

}

void random_spin_state(double (&theta)[N], double (&phi)[N]){
for (int i = 0; i < N; i++){
                double ang1 =  2.0 * pi * random_uniform();
                double ang2 =  acos((2.0 * random_uniform())-1);
                theta[i] = ang1; // This corrects for overcounting at the poles of the sphere. See: http://mathworld.wolfram.com/SpherePointPicking.html
                phi[i] = ang2;
    }
}

void uniform_spin_state(double (&theta)[N], double (&phi)[N]){
    for (int i = 0; i < N; i++) {
        theta[i] = 0.0;
        phi[i] = pi / 2;
    }
}

void magnetization(const double theta[N], const double phi[N], double &Mx, double &My, double &Mz, double &M){
    Mx = 0.0;
    My = 0.0;
    Mz = 0.0;
    for (int i = 0; i < N; i++){
                Mx += cos(theta[i]) * sin(phi[i]);
                My += sin(theta[i]) * sin(phi[i]);
                Mz += cos(phi[i]);
    }
    /*Mx = (Mx/N);
    My = (My/N);
    Mz = (Mz/N);*/
    M = sqrt((Mx*Mx) + (My*My) + (Mz*Mz));
}

void single_spin_M(const double &theta, const double &phi, double &Mx, double &My, double &Mz){
    Mx = (cos(theta) * sin(phi));
    //Mx = Mx/N;
    My = (sin(theta) * sin(phi));
    //My = My/N;
    Mz = cos(phi);
    //Mz = Mz/N;
}

void susceptibility(double &chi, const double M_sum, const double M_sumsquares, const double T){
    chi = (M_sumsquares-(M_sum*M_sum))/(T*T);
}

double single_spin_energy(const int i, const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N]){
    auto e = 0.0;
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[up[i]]) * sin(phi[up[i]]);
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[down[i]]) * sin(phi[down[i]]);
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[left[i]]) * sin(phi[left[i]]);
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[right[i]]) * sin(phi[right[i]]);
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[forwards[i]]) * sin(phi[forwards[i]]);
    e += cos(theta[i]) * sin(phi[i]) * cos(theta[backwards[i]]) * sin(phi[backwards[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[up[i]]) * sin(phi[up[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[down[i]]) * sin(phi[down[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[left[i]]) * sin(phi[left[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[right[i]]) * sin(phi[right[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[forwards[i]]) * sin(phi[forwards[i]]);
    e += sin(theta[i]) * sin(phi[i]) * sin(theta[backwards[i]]) * sin(phi[backwards[i]]);
    e += cos(phi[i]) * cos(phi[up[i]]);
    e += cos(phi[i]) * cos(phi[down[i]]);
    e += cos(phi[i]) * cos(phi[left[i]]);
    e += cos(phi[i]) * cos(phi[right[i]]);
    e += cos(phi[i]) * cos(phi[forwards[i]]);
    e += cos(phi[i]) * cos(phi[backwards[i]]);
    return e * (-J/2.0);
}

double local_energy(const int i, const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N]){
    auto e = 0.0;
    e += single_spin_energy(i, theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(up[i], theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(down[i], theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(left[i], theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(right[i], theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(backwards[i], theta, phi, up, down, left, right, backwards, forwards);
    e += single_spin_energy(forwards[i], theta, phi, up, down, left, right, backwards, forwards);
    return e;
}

void total_energy(const double (&theta)[N], const double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], double &E) {
    E = 0.0;
    for (int i = 0; i < N; i++){
                E += single_spin_energy(i, theta, phi, up, down, left, right, backwards, forwards);
    }
    //E = E/N;
}

void heat_cap(double &C, const double E_sum, const double E_sumsquares, const double T){
    C = (E_sumsquares - (E_sum*E_sum))/(T*T);
}

void Binder(const double &M_sumsquares, const double &M_sumfour){

    double U = 1.0 - (M_sumfour / (3 * M_sumsquares * M_sumsquares));

    myfile.open("binder.txt", std::ios::app);
    myfile << U << std::endl;
    myfile.close();
}

void MC_parallel(double (&theta)[N], double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], const double T, double &E, double &E_sum, double &E_sumsquares, double &Mx, double &My, double &Mz, double &M, double &M_sum, double &M_sumsquares, double &M_sumfour, const int b, double &moves_accepted) {
    E_sum = 0.0;
    E_sumsquares = 0.0;
    M_sum = 0.0;
    M_sumsquares = 0.0;
    M_sumfour = 0.0;
    moves_accepted = 0.0;
    auto counter = 0;
    double E1, E2, Mx1, Mx2, My1, My2, Mz1, Mz2, theta_trial[N] = {}, phi_trial[N] = {};
    double dE;
    for(int i = 0; i < N; i++){
        theta_trial[i] = theta[i];
        phi_trial[i] = phi[i];
    }
    for (int a = 0; a < sweeps; a++) {
        for (int i = 0; i < N; i++) {

                    //total_energy(theta, phi, up, down, left, right, backwards, forwards, E1);
                    E1 = local_energy(i, theta, phi, up, down, left, right, backwards, forwards);

                    theta_trial[i] = 2.0 * pi * random_uniform();
                    phi_trial[i] = std::acos((2.0 * random_uniform()) - 1);


                    E2 = local_energy(i, theta_trial, phi_trial, up, down, left, right, backwards, forwards);
                    //total_energy(theta_trial, phi_trial, up, down, left, right, backwards, forwards, E2);


                    dE = E2 - E1;

                    if (random_uniform() <= exp(std::min(0.0, -dE / T))) {
                        E = E + dE;
                        single_spin_M(theta[i], phi[i], Mx1, My1, Mz1);
                        single_spin_M(theta_trial[i], phi_trial[i], Mx2, My2, Mz2);
                        Mx += (Mx2 - Mx1);
                        My += (My2 - My1);
                        Mz += (Mz2 - Mz1);
                        M = sqrt((Mx*Mx) + (My*My) + (Mz*Mz));
                        theta[i] = theta_trial[i];
                        phi[i] = phi_trial[i];
                        moves_accepted += 1.0;
                    }
                    else{
                        theta_trial[i] = theta[i];
                        phi_trial[i] = phi[i];
                    }


                    if (a % 5 == 0) {
                        E_sum += E;
                        E_sumsquares += (E * E);
                        M_sum += M;
                        M_sumsquares += (M * M);
                        M_sumfour += (M * M * M * M);
                        counter++;
                    }
        }
    }
    E_sum = E_sum / counter;
    E_sumsquares = E_sumsquares / counter;
    M_sum = M_sum / counter;
    M_sumsquares = M_sumsquares / counter;
    M_sumfour = M_sumfour / counter;
}

void thermalize_parallel(double (&theta)[N], double (&phi)[N], const int (&up)[N], const int (&down)[N], const int (&left)[N], const int (&right)[N], const int (&backwards)[N], const int (&forwards)[N], const double T){
    double E1, E2, dE, theta_trial[N] = {}, phi_trial[N] = {};
    for(int i = 0; i < N; i++){
        theta_trial[i] = theta[i];
        phi_trial[i] = phi [i];
    }
    for (int a = 0; a < thermalsteps; a++) {
        for (int i = 0; i < N; i++) {

                    E1 = local_energy(i, theta, phi, up, down, left, right, backwards, forwards);

                    theta_trial[i] = 2.0 * pi * random_uniform();
                    phi_trial[i] = std::acos((2.0 * random_uniform()) - 1);



                    E2 = local_energy(i, theta_trial, phi_trial, up, down, left, right, backwards, forwards);

                    dE = E2 - E1;
                    if (random_uniform() <= exp(std::min(0.0, -dE / T))) {
                        theta[i] = theta_trial[i];
                        phi[i] = phi_trial[i];
                    }
                    else{
                        theta_trial[i] = theta[i];
                        phi_trial[i] = phi[i];
                    }
        }
    }
}

void iterate_T(const double &T){
    myfile.open("temp.txt", std::ios::app);
    myfile << T << std::endl;
    myfile.close();
    //T += T_step;
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
