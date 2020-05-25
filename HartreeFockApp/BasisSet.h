#pragma once
#include <vector>
#include <iostream>
#include <armadillo>
#include <chrono>
#include <future>
#include <fstream>
#include "Nucleons.h"
#include "ContractedGTO.h"

class BasisSet
{
public:


public:

	BasisSet() {};

	void createBasisSet(const std::vector<Nucleon> & nucleons);
	void createTestBasisSet(const std::vector<Nucleon> & nucleons, double l, double s);

	void calculateIntegrals();

	std::vector<std::vector<std::tuple<double, double, Position>>> loadAOs(int n);
	void saveAOs(const std::vector<std::vector<std::tuple<double, double, Position>>>&, int n);

	double getSize();

	double getH(int i, int j);
	double getS(int i, int j);
	double getDI(int i, int j, int k, int l);

	double getFunctionValue(int n, const Position& p);

private:

	std::vector<ContractedGTO> _basisFucntion;

	std::vector<std::vector<double>> _S_Matrix;

	std::vector<std::vector<double>> _H_Matrix;

	std::vector<std::vector<std::vector<std::vector<double>>>> _D_Matrix;

	std::vector<Nucleon> _nucleons;

private:

	void calculateOneElectronHamiltonians();
	void calculateOverlapIntegrals();
	void calculateTwoElectronIntegrals();

};