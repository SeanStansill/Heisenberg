//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <fstream>


const int L = 32;
const int thermalsteps = 2000000; //2000000
const int nsteps = 50000;
const double J = 1.0;
const double k_B = 1.0;
const double pi = 3.141592654;
const double T_init = 0.0;
const double T_step = 0.05;
const int T_intervals = 81;
const int N = L*L*L;
const int latticesweeps = 500000;
const double exp1 = 2.718281828;
const int steps = 10000;