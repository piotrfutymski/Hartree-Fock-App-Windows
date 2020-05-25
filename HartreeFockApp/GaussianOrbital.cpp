#include "GaussianOrbital.h"

void GaussianOrbital::setalfa(double a)
{
	alfa = a;
	this->recalculateNormalizationParam();
}

Position GaussianOrbital::getR()
{
	return R0;
}

double GaussianOrbital::getalfa()
{
	return alfa;
}

double GaussianOrbital::getNormalizationParam()
{
	return normalizationParam;
}

double GaussianOrbital::calculateOverlapIntegral(const GaussianOrbital & right)const
{
	double b = this->alfa + right.alfa;
	double A = 4.0 * (this->alfa * right.alfa) / (b * b);
	auto R = this->R0 - right.R0;

	return pow(A, 0.75) * exp(-(this->alfa * right.alfa) * R.r2 / b);
		
}

double GaussianOrbital::calculateHIntegral(const GaussianOrbital & right, const std::vector<Nucleon>& nucleons)const
{
	return this->calculateKineticIntegral(right) - this->calculateNucleonIntegral(right, nucleons);
}

double GaussianOrbital::calculateHIntegral(const GaussianOrbital& right, const std::vector<Nucleon>& nucleons, double S)const
{
	return this->calculateKineticIntegral(right, S) - this->calculateNucleonIntegral(right, nucleons, S);
}

double GaussianOrbital::calculateKineticIntegral(const GaussianOrbital & right)const
{
	return this->calculateKineticIntegral(right, this->calculateOverlapIntegral(right));
}

double GaussianOrbital::calculateKineticIntegral(const GaussianOrbital& right, double S)const
{
	double b = this->alfa + right.alfa;
	double C = this->alfa * right.alfa / b;
	auto R = this->R0 - right.R0;
	return(C * (3.0 - 2 * C * R.r2) * S);
}

double GaussianOrbital::calulateTwoElectronIntegral(const GaussianOrbital & a, const GaussianOrbital & b, const GaussianOrbital & c, const GaussianOrbital & d)
{
	return GaussianOrbital::calulateTwoElectronIntegral(a, b, c, d, a.calculateOverlapIntegral(c), b.calculateOverlapIntegral(d));
}

double GaussianOrbital::calulateTwoElectronIntegral(const GaussianOrbital& a, const GaussianOrbital& b, const GaussianOrbital& c, const GaussianOrbital& d, double Spq, double Srs)
{
	double A = a.alfa + c.alfa;
	double B = b.alfa + d.alfa;
	double AB = A + B;

	auto Rk = (1.0/A)*(a.alfa * a.R0 + c.alfa * c.R0);
	auto Rl = (1.0/B)*(b.alfa * b.R0 + d.alfa * d.R0);

	auto R = Rk - Rl;

	return 2.0 / sqrt(PI) * (sqrt(A) * sqrt(B)) / sqrt(AB) * BoysCalculator::boys(A * B / AB * R.r2) * Spq * Srs;
}

double GaussianOrbital::calculateNucleonIntegral(const GaussianOrbital & right, const std::vector<Nucleon>& nucleons)const
{
	return this->calculateNucleonIntegral(right, nucleons, this->calculateOverlapIntegral(right));
}

double GaussianOrbital::calculateNucleonIntegral(const GaussianOrbital& right, const std::vector<Nucleon>& nucleons, double S)const
{
	double res = 0;
	for (size_t i = 0; i < nucleons.size(); i++)
	{
		res += calculateCulombIntegral(right, nucleons[i].p, S) * nucleons[i].charge;
	}
	return res;
}



void GaussianOrbital::recalculateNormalizationParam()
{
	normalizationParam = pow(2 * alfa / PI, 0.75);
}


double GaussianOrbital::calculateCulombIntegral(const GaussianOrbital& right, const Position& p) const
{
	return this->calculateCulombIntegral(right, p, this->calculateOverlapIntegral(right));
}

double GaussianOrbital::calculateCulombIntegral(const GaussianOrbital& right, const Position& p, double S) const
{
	double A = this->alfa + right.alfa;
	auto Rk = (1.0 / A) * (this->alfa * this->R0 + right.alfa * right.R0);
	auto R = p - Rk;
	return 2.0 * sqrt(A / PI) * BoysCalculator::boys(A * R.r2) * S;
}


std::complex<double> GaussianOrbital::F_X(const Position & p)
{
	auto pos = p - R0;
	return std::complex<double>(this->getNormalizationParam()*exp(-alfa * pos.r2), 0);
}
