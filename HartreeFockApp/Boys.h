#pragma once
#include <cmath>
#include "Constatns.h"
#include <map>
#include <vector>

class BoysCalculator
{
static std::map<double, double> boysTab;

public:

	static int silnia(int n);
	static double boysF1(double x);
	static double boys(double x);
	static double boysHighz(int m, double z, int N);
	static double boysGaussJacobi(int m, double z);
};















