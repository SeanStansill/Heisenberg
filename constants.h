//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//
#include <iostream>
#include <fstream>


const int L = 16;
const int thermalsteps = 40000; //2000000
const int nsteps = 200000;
const double J = 1.0;
const double k_B = 1.0;
const double pi = 3.141592654;
const double T_init = 0.05;
const double T_step = 0.05;
const int T_intervals = 60;
const int N = L*L*L;
const int steps = 10000;
const int sweeps = 1000;
const int nthreads = 8;
const int parallel_sweeps = sweeps/nthreads;