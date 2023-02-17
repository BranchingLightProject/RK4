#ifndef __BRANCHED_FLOW_H
#define __BRANCHED_FLOW_H

#include<iostream>
#include<cmath>
#include<complex>
#include<string>

#include"inter/stdafx.h"
#include"inter/interpolation.h"

#include"constants.h"
#include"random64.h"
#include"csv_handler.h"

typedef std::complex<double> c_double;
typedef unsigned long long u_l_long;

class BranchedFlow{
    private:
        CRandom *ran = NULL;
        alglib::spline2dinterpolant inter;

        double *film = NULL;
        double *scint = NULL;

        c_double phi[Ny] = {0.0};
        std::vector<double> corr_funct;
        std::vector<std::vector<double>> corr_funct_2D;

        double dr = 0.0, corr_length = 0.0;
        double dx = 0.0, dy = 0.0;
        int n_bins = 0, nx_bins = 0, ny_bins = 0;
    public:
        #if POTENTIAL_SOURCE == 0
        BranchedFlow(std::string filename, u_l_long seed, double shift=0.0);
        #else
        BranchedFlow(unsigned long long seed);
        #endif

        ~BranchedFlow();

        void initialize(int initial=0);

        void rk4_solve(void);

        void corr_solve(int bins=500);
        void corr_solve_2D(int bins_x=100, int bins_y=100);
        double potential_strength(void);
        double potential_average(void);

        void save_film(std::string filename);
        void save_potential(std::string filename, bool effective=true);
        void save_scint(std::string filename);
        void save_corr(std::string filename);
        void save_corr_2D(std::string filename);

        double correlation_length(void){return corr_length;}
};

#endif // __BRANCHED_FLOW_H
