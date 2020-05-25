#pragma once

#include <complex>
#include <cmath>
#include <vector>
#include <functional>
#include "Position.h"

class Orbital
{
public:

	virtual std::complex<double> F_X(const Position &) = 0;

private:

public:

	static std::function<std::complex<double>(const Position &)> getGaussianFunctionR(int n, int l, const Position& R0, double b);

	static std::function<std::complex<double>(const Position &)> getTrueFunctionR(int n, int l, const Position& R0);

	static std::function<std::complex<double>(const Position &)> getSphericalFunction(int l, int m, const Position& R0);

	static double leguerrePolynomial(double x, int k, int p);

	static std::vector<double> leguarrePolynomialCo(int k, int p);

	static double get_x_n_e_x_2_Calc(double beta, int n);
};