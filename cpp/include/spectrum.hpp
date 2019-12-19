#ifndef _SPECTRUM_HPP
#define _SPECTRUM_HPP



#include <cmath>
#include <limits>
#include <vector>

#include "gsl_const_cgsm.h"
#include "util.hpp"


namespace Spectrum {
constexpr const double double_h_over_c2 = 2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H / (GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT);
constexpr const double h_over_kB = GSL_CONST_CGSM_PLANCKS_CONSTANT_H / GSL_CONST_CGSM_BOLTZMANN;
constexpr const double double_h_c2 = 2. * GSL_CONST_CGSM_PLANCKS_CONSTANT_H * GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT;
constexpr const double ch_over_kB = GSL_CONST_CGSM_SPEED_OF_LIGHT * h_over_kB;
double Planck_nu(double T, double nu);
double Planck_lambda(double T, double lambda);

double Planck_nu1_nu2(double T, double nu1, double nu2, double tol=std::sqrt(std::numeric_limits<double>::epsilon()));

class DiskGR {
private:
    const double kerr;
    const double x0;
    const double x1;
    const double x2;
    const double x3;
    const double a0;
    const double a1;
    const double b0;
    const double b1;
    const double c0;
    const double c1;
    inline double a(const double x) { return a0 * (std::log(x-x1) - a1); }
    inline double b(const double x) { return b0 * (std::log(x-x2) - b1); }
    inline double c(const double x) { return c0 * (std::log(x-x3) - c1); }
public:
    DiskGR(double kerr);
    double T(double r, double Mx, double Mdot) const;
};

double T_GR(double r1, double ak, double Mx, double Mdot);
} // namespace Spectrum


#endif // _SPECTRUM_HPP
