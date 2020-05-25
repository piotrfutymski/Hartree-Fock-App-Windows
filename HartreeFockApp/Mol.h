#pragma once
#include <armadillo>
#include "Nucleons.h"
#include "BasisSet.h"
#include "EigenSolver.h"
#include <chrono>
#include <future>

class Mol
{
public:

	Mol(const std::vector<Nucleon>& n);
	Mol();
	Mol(int charge);

	void moveNucleon(int n, const Position& delta);
	void setMOcount(int c);


	void initBasisSet();
	void calculateIntegrals();
	void HFProcedure(int n = 0);

	void HF_TO_Divergance(double delta);
	void HF_TO_ElapsedTime(double t);

	double getElectronicEnergy();
	double getMoleculeEnergy();

	std::vector<double> getMolecularCoeficents(int m);
	double countMolecularFunction(int m, const Position & p);
	double countMolecularFunction(int m, double x);

private:

	std::vector<Nucleon> _nucleons;
	int _MOcount;
	BasisSet _basisSet;
	bool _integralsCalculated{ false };
	bool _setInitiated{ false };

	arma::mat _FMatrix;
	arma::mat _SMatrix;
	arma::mat _PMatrix;

	arma::mat _orbitalCoeficients[5];

	double _electronicEnergy;

private:

	void recalculateP();
	void recalculateF();
	void recalculateEnergy();

	arma::vec getCoeficientIntegrals(const arma::mat& left, const arma::mat& right);

	void PullyMixing(int n);

};