//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <fstream>


const int L = 16;
const int thermalsteps = 200000;
const int nsteps = 10000;
const double J = 1.0;
const double k_B = 1.0;
const double pi = 3.141592654;
const double T_init = 0.1;
const double T_step = 0.1;
const int T_intervals = 20;
const int N = L*L*L;
