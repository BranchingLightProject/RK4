#include "branched_flow.h"

#if POTENTIAL_SOURCE == 0
BranchedFlow::BranchedFlow(std::string filename, u_l_long seed, double shift) {
    ran = new CRandom(seed);
    scint = new double[Nx];
    film = new double[Sx * Sy];

    // Potential
    CsvHandler file(filename);

    int Lx = file.get_rows();
    int Ly = file.get_cols();

    double V[Lx * Ly];

    double x[Lx], y[Ly];

    double x_factor = Lx_real / (double)Lx, y_factor = Ly_real / (double)Ly;

    for (int ix = 0; ix < Lx; ix++) {
        x[ix] = (double)ix * x_factor;
        for (int iy = 0; iy < Ly; iy++) {
            y[iy] = (double)iy * y_factor;
            V[Lx * iy + ix] = file.content[ix][Ly - 1 - iy] - shift;
        }
    }

    alglib::real_1d_array a_x, a_y, a_V;

    a_x.setcontent(Lx, x);
    a_y.setcontent(Ly, y);
    a_V.setcontent(Lx * Ly, V);

    alglib::spline2dbuildbicubicv(a_x, Lx, a_y, Ly, a_V, 1, inter);
}
#else
BranchedFlow::BranchedFlow(unsigned long long seed) {
    ran = new CRandom(seed);
    scint = new double[Nx];
    film = new double[Sx * Sy];
}
#endif

BranchedFlow::~BranchedFlow() {
    delete ran;
    delete[] scint;
    delete[] film;
}

void BranchedFlow::initialize(int initial) {
    if (initial == 0) {
        for (int j = 0; j < Ny; j++)
            phi[j] = std::exp(-0.5 * (l_d_eta * j - beam_mu) * (l_d_eta * j - beam_mu) / (beam_sigma * beam_sigma));
    } else if (initial == 1) {
        for (int j = 0; j < Ny; j++)
            phi[j] = M_1_PI * (std::atan(25.0 * (l_d_eta * j - 0.2)) - std::atan(25.0 * (l_d_eta * j - 4.8)));
    }
}

#if POTENTIAL_SOURCE == 0
// Potential interpolated from file
#define potential(x, y) (alglib::spline2dcalc(inter, x, y))

#else
double period = 1.0;  // mm

double factor_x = 2.0 * M_PI / period;
double factor_y = 2.0 * M_PI / period;
/**
 * Potential as a function of x and y. The general form is
 * K0^2*( (A*f(x,y) + N0^2) - N0^2) = K0^2*A*f(x, y)
 *
 * Please note that x and y are NOT dimensionless. They must be in millimeters.
 *
 * Used potentials:
 *  K2*0.52*M_2_PI*std::atan( (2.0*((x)-6.5) + (y)) ); // angled glass
 *  K2*0.1*(std::sin(factor_x*(x)) + std::sin(factor_y*(y))); // simple periodic potential
 *  K2*0.1*(std::sin(factor_x*(x)-factor_y*(y)) + std::sin(factor_y*(y))); // angled periodic potential
 *  K2*0.1*(0.75 + 0.5*(x)/Lx_real)*(std::sin(factor_x*(x)-5*factor_y*(y)) + std::sin(factor_y*((y)-100))); // angled and twisted p.p.
 */
#define potential(x, y) (0.52 * M_2_PI * std::atan((((x)-5.0) + (y))))  // angled glass)
#endif

// Differential equation to be solved numerically
#define paraxial_equation(uy_1, uy, uy1, v) (C_eq * ((uy1 - 2.0 * uy + uy_1) * o_d_eta2 + C_v * v * uy))

// Evolve the initial profile
void BranchedFlow::rk4_solve(void) {
    static const double d_Ny = (double)Ny, Px = 100.0 / (double)Nx;
    double p = 0.0;  // progress

    static const int max_boundary = 1500;
    static const double factor = 1.0 / (double)max_boundary;
    double f = 0.0;

    static c_double k[4][Ny] = {0.0};

    for (int i = 0, film_i = 0; i < Nx; i++, film_i--, p += Px) {
        /* ~~ k1 ~~ */
        for (int j = 1; j < Ny - 1; j++) {
            c_double uy_1 = phi[j - 1];
            c_double uy = phi[j];
            c_double uy1 = phi[j + 1];
            double v = potential(l_d_xi * i, l_d_eta * j);

            k[0][j] = d_xi * paraxial_equation(uy_1, uy, uy1, v);
        }

        /* ~~ k2 ~~ */
        for (int j = 1; j < Ny - 1; j++) {
            c_double uy_1 = phi[j - 1] + (k[0][j - 1] / 2.0);
            c_double uy = phi[j] + (k[0][j] / 2.0);
            c_double uy1 = phi[j + 1] + (k[0][j + 1] / 2.0);
            double v = potential(l_d_xi * (i + 0.5), l_d_eta * (j + 0.5));

            k[1][j] = d_xi * paraxial_equation(uy_1, uy, uy1, v);
        }

        /* ~~ k3 ~~ */
        for (int j = 1; j < Ny - 1; j++) {
            c_double uy_1 = phi[j - 1] + (k[1][j - 1] / 2.0);
            c_double uy = phi[j] + (k[1][j] / 2.0);
            c_double uy1 = phi[j + 1] + (k[1][j + 1] / 2.0);
            double v = potential(l_d_xi * (i + 0.5), l_d_eta * (j + 0.5));

            k[2][j] = d_xi * paraxial_equation(uy_1, uy, uy1, v);
        }

        /* ~~ k4 ~~ */
        for (int j = 1; j < Ny - 1; j++) {
            c_double uy_1 = phi[j - 1] + k[2][j - 1];
            c_double uy = phi[j] + k[2][j];
            c_double uy1 = phi[j + 1] + k[2][j + 1];
            double v = potential(l_d_xi * (i + 1.0), l_d_eta * (j + 1.0));

            k[3][j] = d_xi * paraxial_equation(uy_1, uy, uy1, v);
        }

        double I_prom = 0.0, I2_prom = 0.0;

        // Evolve
        for (int j = 0, film_j = 0; j < Ny; j++, film_j--) {
            double norm = std::norm(phi[j]);

            if (!(film_i + film_j)) {
                film[(i / (x_save_step + 1)) * Sy + j / (y_save_step + 1)] = norm;
                film_j = y_save_step;
            }

            I_prom += norm;
            I2_prom += norm * norm;

            phi[j] = phi[j] + (k[0][j] + 2.0 * k[1][j] + 2.0 * k[2][j] + k[3][j]) / 6.0;
        }

        // Boundary conditions
        phi[0] = 0.0;
        phi[Ny - 1] = 0.0;
        f = factor;
        for (int b = 1; b < max_boundary; b++, f += factor) {
            phi[b] *= f;
            phi[Ny - 1 - b] *= f;
        }

        // Scintillation
        scint[i] = -1.0 + d_Ny * I2_prom / (I_prom * I_prom);

        if (film_i == 0) {
            film_i = x_save_step;

            // Progress bar
            printf("\r%.2f %%", p);
            fflush(stdout);
        }
    }
}
#undef paraxial_equation

// Find c(dr)
void BranchedFlow::corr_solve(int bins) {
    double products[bins];
    products[0] = 1.0;
    double v_prom = 0.0, v2_prom = 0.0;
    double Px = 100.0 / (double)bins, p = 0.0;  // progress

    dr = 4.0 / bins;  // mm/bins
    n_bins = bins;

    // Average product by Montecarlo sampling
    for (int i = 1; i < bins; i++, p += Px) {
        double R = i * dr;

        for (int t = 0; t < CORR_STEPS; t++) {
            double x1 = Lx_real * ran->r(), y1 = Ly_real * ran->r();

            double phi = 2.0 * M_PI * ran->r();

            double x2 = std::fmod(x1 + R * std::cos(phi) + Lx_real, Lx_real);
            double y2 = std::fmod(y1 + R * std::sin(phi) + Ly_real, Ly_real);

            products[i] += potential(x1, y1) * potential(x2, y2);
        }

        // Progress bar
        printf("\r%.2f %%", p);
        fflush(stdout);
    }

    // Average potential and potential squared
    for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++) {
            double v = potential(lambda * i * d_xi * x_save_step, lambda * j * d_eta * y_save_step);
            v_prom += v;
            v2_prom += v * v;
        }

    v_prom /= (double)(Sx * Sy);
    v2_prom /= (double)(Sx * Sy);

    double v_prom2 = v_prom * v_prom;

    corr_funct.clear();
    corr_funct.push_back(1.0);
    double c_l = dr;

    // Load to vector and calculate correlation length
    for (int i = 1; i < bins; i++) {
        double average = products[i] / (double)CORR_STEPS;
        double corr = (average - v_prom2) / (v2_prom - v_prom2);

        corr_funct.push_back(corr);
        c_l += corr * dr;
    }

    corr_length = c_l;
}

void BranchedFlow::corr_solve_2D(int bins_x, int bins_y) {
    double products[bins_x][bins_y];
    double v_prom = 0.0, v2_prom = 0.0;
    double Px = 100.0 / (double)bins_x, p = 0.0;  // progress

    nx_bins = bins_x;
    ny_bins = bins_y;

    dx = 2.0 / bins_x;  // mm/bins
    dy = 2.0 / bins_y;  // mm/bins

    // Average product by Montecarlo sampling
    for (int ix = 0; ix < bins_x; ix++, p += Px) {
        double Rx = ix * dx;

        for (int iy = 0; iy < bins_y; iy++) {
            double Ry = iy * dy;

            for (int t = 0; t < CORR_STEPS; t++) {
                double x1 = Lx_real * ran->r(), y1 = Ly_real * ran->r();

                double phi = 2.0 * M_PI * ran->r();

                double x2 = std::fmod(x1 + Rx * std::cos(phi) + Lx_real, Lx_real);
                double y2 = std::fmod(y1 + Ry * std::sin(phi) + Ly_real, Ly_real);

                products[ix][iy] += potential(x1, y1) * potential(x2, y2);
            }
        }

        // Progress bar
        printf("\r%.2f %%", p);
        fflush(stdout);
    }

    // Average potential and potential squared
    for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++) {
            double v = potential(lambda * i * d_xi * x_save_step, lambda * j * d_eta * y_save_step);
            v_prom += v;
            v2_prom += v * v;
        }

    v_prom /= (double)(Sx * Sy);
    v2_prom /= (double)(Sx * Sy);

    double v_prom2 = v_prom * v_prom;

    corr_funct_2D.clear();
    double c_l = 0.0;

    // Load to vector and calculate correlation length
    for (int ix = 0; ix < bins_x; ix++) {
        std::vector<double> line;
        for (int iy = 0; iy < bins_y; iy++) {
            double average = products[ix][iy] / CORR_STEPS;
            double corr = (average - v_prom2) / (v2_prom - v_prom2);

            line.push_back(corr);
        }
        corr_funct_2D.push_back(line);
    }
}

// sqrt( average(V^2) )/2E
double BranchedFlow::potential_strength(void) {
    double v0 = 0.0;

    for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++) {
            double v = potential(l_d_xi * i * x_save_step, l_d_eta * j * y_save_step);
            v0 += v * v;
        }

    v0 /= (double)(Sx * Sy);
    v0 = std::sqrt(v0);

    // TODO: fix units
    v0 /= 2 * K2 * N2 * 1e-8;

    return v0;
}

// average(V)
double BranchedFlow::potential_average(void) {
    double v0 = 0.0;

    for (int i = 0; i < Sx; i++)
        for (int j = 0; j < Sy; j++) {
            v0 += potential(l_d_xi * i * x_save_step, l_d_eta * j * y_save_step);
        }

    v0 /= (double)(Sx * Sy);

    return v0;
}

/**
 * Save used potential
 *
 * @param filename: name of file
 * @param effective: save potential as experienced by the system. (This multiplies by the constant of the equation).
 */
void BranchedFlow::save_potential(std::string filename, bool effective) {
    std::ofstream file(filename);

    double factor = 1.0;

    if (effective) factor = C_v;

    for (int i = 0; i < Sx; i++) {
        for (int j = 0; j < Sy; j++)
            file << lambda * i * d_xi * x_save_step << ',' << lambda * j * d_eta * y_save_step << ','
                 << factor * potential(lambda * i * d_xi * x_save_step, lambda * j * d_eta * y_save_step) << '\n';
        file << '\n';
    }

    file << std::endl;
    file.close();
}

#undef potential

void BranchedFlow::save_film(std::string filename) {
    std::ofstream file(filename);

    double x_mult = (double)Nx / (double)Sx, y_mult = (double)Ny / (double)Sy;

    for (int i = 0; i < Sx; i++) {
        for (int j = 0; j < Sy; j++)
            file << lambda * i * d_xi * x_mult << ',' << lambda * j * d_eta * y_mult << ','
                 << std::sqrt(film[i * Sy + j]) << '\n';
        file << '\n';
    }

    file << std::endl;
    file.close();
}

void BranchedFlow::save_scint(std::string filename) {
    std::ofstream file(filename);

    for (int i = 0; i < Nx; i++)
        file << lambda * i * d_xi << ',' << scint[i] << '\n';

    file.close();
}

void BranchedFlow::save_corr(std::string filename) {
    std::ofstream file(filename);

    for (int i = 0; i < n_bins; i++)
        file << i * dr << ',' << corr_funct[i] << '\n';

    file.close();
}

void BranchedFlow::save_corr_2D(std::string filename) {
    std::ofstream file(filename);

    for (int i = nx_bins - 1; i >= 0; i--) {
        for (int j = ny_bins - 1; j >= 0; j--)
            file << -1.0 * i * dx << ',' << -1.0 * j * dy << ',' << corr_funct_2D[i][j] << '\n';
        for (int j = 0; j < ny_bins; j++)
            file << -1.0 * i * dx << ',' << j * dy << ',' << corr_funct_2D[i][j] << '\n';
        file << '\n';
    }

    for (int i = 0; i < nx_bins; i++) {
        for (int j = ny_bins - 1; j >= 0; j--)
            file << i * dx << ',' << -1.0 * j * dy << ',' << corr_funct_2D[i][j] << '\n';
        for (int j = 0; j < ny_bins; j++)
            file << i * dx << ',' << j * dy << ',' << corr_funct_2D[i][j] << '\n';
        file << '\n';
    }

    file.close();
}
