#include "Orbital.h"

std::function<std::complex<double>(const Position&)> Orbital::getTrueFunctionR(int n, int l, const Position & R0)
{
	return [=](const Position& r) {
		auto ra = r - R0;
		auto rn = ra.r / (double)n;
		return std::complex<double>(pow(rn, l)*exp(-rn)*leguerrePolynomial(2 * rn, 2l + 1, n - l - 1));
	};
}

std::function<std::complex<double>(const Position&)> Orbital::getSphericalFunction(int l, int m, const Position & R0)
{
	return std::function<std::complex<double>(const Position&)>();
}

double Orbital::leguerrePolynomial(double x, int k, int p)
{
	if (p > 3 || p < 0)
		throw std::exception("QUANTUM NUMBERS - error with quantum number (max allowed n is 4)");
	else if (p == 0)
		return 1;
	else if (p == 1)
		return -x + k + 1;
	else if (p == 2)
		return x * x / 2 - (k + 2)*x + (k + 1)*(k + 2) / 2;
	else
		return -x * x*x / 6 + (k + 3)*x*x / 2 - (k + 2)*(k + 3)*x / 2 + (k + 1)*(k + 2)*(k + 3) / 6;
}

std::vector<double> Orbital::leguarrePolynomialCo(int k, int p)
{
	if (p > 3 || p < 0)
		throw std::exception("QUANTUM NUMBERS - error with quantum number (max allowed n is 4)");
	else if (p == 0)
		return { 1.0 };
	else if (p == 1)
		return { 1.0 + k, -2 };
	else if (p == 2)
		return { (k + 1.0)*(k + 2.0) / 2.0 , -2*(k + 2.0) , 2.0 };
	else
		return { (k + 1.0)*(k + 2.0)*(k + 3.0) / 6.0 ,  (k + 2.0)*(k + 3.0), (k + 3.0) *2.0, -8.0 / 6.0 };
}

double Orbital::get_x_n_e_x_2_Calc(double beta, int n)
{
	if (n == 0)
		return sqrt(PI / beta);
	else if (n == 1)
		return 1.0 / ( beta);
	else if (n == 2)
		return sqrt(PI / beta) / (2 * beta);
	else if (n == 3)
		return 1.0 / ( beta * beta);
	else
	{
		if (n % 2 == 0)
			return (n - 1)*get_x_n_e_x_2_Calc(beta, n - 2) / (2 * beta);
		else
			return ((n - 1) / 2)*get_x_n_e_x_2_Calc(beta, n - 2) / beta;
	}
}

std::function<std::complex<double>(const Position&)> Orbital::getGaussianFunctionR(int n, int l, const Position& R0, double b)
{
	return [=](const Position& r) {
		auto ra = r - R0;
		auto rn = ra.r / (double)n;
		return std::complex<double>(pow(rn, l)*exp(-b * ra.r * ra.r)*leguerrePolynomial(2 * rn, 2*l + 1, n - l - 1));
	};
}
