//
// Created by Sean Stansill [ll14s26s] on 25/09/2019.
//

const int N = 32;
const int nsteps = 1000;
const double J = 1.21e-21;  //value for iron so Tc ~ 1000K for reference from https://www.southampton.ac.uk/~rpb/thesis/node18.html
const double k_B = 1.38e-23;
const double pi = 3.141592654;
const double T_init = 10;
const double T_step = 10;
const int T_intervals = 200;
const double L = N*N*N;
