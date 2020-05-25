#pragma once
#include "GaussianOrbital.h"

class ContractedGTO : public Orbital
{
public:

	void addPrimitive(double, const GaussianOrbital&);

	double calculateOverlapIntegral(const ContractedGTO& right) const;
	double calculateHIntegral(const ContractedGTO& right, const std::vector<Nucleon>& nucleons)const;
	double calculateHIntegral(const ContractedGTO& right, const std::vector<Nucleon>& nucleons, double S)const;
	static double calulateTwoElectronIntegral(const ContractedGTO& a, const ContractedGTO& b, const ContractedGTO& c, const ContractedGTO& d);
	static double calulateTwoElectronIntegral(const ContractedGTO& a, const ContractedGTO& b, const ContractedGTO& c, const ContractedGTO& d, double Spq, double Srs);

private:

	std::vector<std::pair<double, GaussianOrbital>> _primitives;

	// Odziedziczono za poœrednictwem elementu Orbital
public:
	virtual std::complex<double> F_X(const Position&) override;

};