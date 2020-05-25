#pragma once
#include "Orbital.h"
#include "Nucleons.h"
#include "Boys.h"
#include <map>

class GaussianOrbital : public Orbital
{
public:

	GaussianOrbital(const Position & R, double a) :R0{ R}, alfa{a}
	{
		this->recalculateNormalizationParam();
	}

	void setalfa(double a);

	Position getR();
	double getalfa();
	double getNormalizationParam();

	// Integrals

	double calculateOverlapIntegral(const GaussianOrbital & right) const;
	//
	double calculateNucleonIntegral(const GaussianOrbital & right, const std::vector<Nucleon> & nucleons)const;
	double calculateNucleonIntegral(const GaussianOrbital & right, const std::vector<Nucleon> & nucleons, double S)const;
	double calculateHIntegral(const GaussianOrbital & right, const std::vector<Nucleon> & nucleons)const;
	double calculateHIntegral(const GaussianOrbital & right, const std::vector<Nucleon> & nucleons, double S)const;
	double calculateKineticIntegral(const GaussianOrbital & right)const;
	double calculateKineticIntegral(const GaussianOrbital & right, double S)const;
	//
	static double calulateTwoElectronIntegral(const GaussianOrbital &a, const GaussianOrbital &b, const GaussianOrbital &c, const GaussianOrbital &d);
	static double calulateTwoElectronIntegral(const GaussianOrbital &a, const GaussianOrbital &b, const GaussianOrbital &c, const GaussianOrbital &d, double Spq, double Srs);

private:

	Position R0;
	double alfa;
	double normalizationParam;

private:

	void recalculateNormalizationParam();

	//Culomb
	double calculateCulombIntegral(const GaussianOrbital& right, const Position& p)const;
	double calculateCulombIntegral(const GaussianOrbital& right, const Position& p, double S)const;


public:
	// Inherited via Orbital
	virtual std::complex<double> F_X(const Position &) override;

};

