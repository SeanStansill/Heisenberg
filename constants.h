//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <fstream>


const int L = 4;
const int thermalsteps = 500000;
const int nsteps = 10000;
const double J = 1.0;
const double k_B = 1.0;
const double pi = 3.141592654;
const double T_init = 0.005;
const double T_step = 0.005;
const int T_intervals = 40;
const int N = L*L*L;

/*const int L = 4; This is the next attempt
const int thermalsteps = 2000000;
const int nsteps = 500000;
const double J = 1.0;
const double k_B = 1.0;
const double pi = 3.141592654;
const double T_init = 0.01;
const double T_step = 0.01;
const int T_intervals = 100;
const int N = L*L*L;*/