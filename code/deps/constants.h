#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include<complex>

typedef std::complex<double> c_double;

/**
 * Determines wether the potential is defined from a function or from an external file.
 * 
 * Options:
 *  * 0: From a file
 *  * 1: From a function 
 */
#define POTENTIAL_SOURCE 0

// Physical constants
static const double N = 1.33; // n_water
static const double N2 = N*N; // n_water^2
static const double lambda = 532.0e-6; // (mm)
static const double K = 2*M_PI/lambda; // wavenumber (1/mm)
static const double K2 = K*K; // wavenumber (1/mm^2)
static const c_double j = {0.0, 1.0};

// Real lengths
static const double Lx_real = 10.0; // mm (= 1 cm)
static const double Ly_real = 5.0; // mm (= 0.5 cm)

// Equation constant
static const double C = 5.983269e-02;

static const c_double C_eq = -1.0*j*C;
static const double C_v = 4.0*M_PI*M_PI;

// Algorithm constants
static const int Nx = 150375; // nodes in x
static const int Ny = 76845; // nodes in y

static const int Sx = 2048; // sampled nodes in x
static const int Sy = 1024; // sampled nodes in y

static const double d_xi = 1.250000e-01; // d_xi (mm)

static const double d_eta = 1.223036e-01; // (mm)
static const double o_d_eta2 = 6.685309e+01; // d_eta^(-2) (mm^-2)

static const double l_d_xi = lambda*d_xi; // lambda*d_xi
static const double l_d_eta = lambda*d_eta; // lambda*d_eta

// Beam & plane wave constants
static const double beam_sigma = 0.125; // Beam half half-width (mm)
static const double beam_mu = (Ly_real/2.0); // Beam center (mm)

// Printing constants
static const int x_save_step = Nx/Sx;
static const int y_save_step = Ny/Sy;

// Correlation
static const int CORR_STEPS = 10000; // number of tests done for c(dr)

#endif // __CONSTANTS_H