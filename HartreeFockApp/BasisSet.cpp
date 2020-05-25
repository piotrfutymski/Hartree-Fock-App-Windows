#include "BasisSet.h"
#include "Mol.h"


void BasisSet::createBasisSet(const std::vector<Nucleon>& nucleons)
{
	_nucleons = nucleons;
	for (int i = 0; i < nucleons.size(); i++)
	{
		auto date = this->loadAOs(_nucleons[i].charge);
		for (int j = 0; j < date.size(); j++)
		{
			_basisFucntion.push_back({});
			for (int k = 0; k < date[j].size(); k++)
			{
				(_basisFucntion.end() - 1)->addPrimitive(
					std::get<1>(date[j][k]), {
						 nucleons[i].p + std::get<2>(date[j][k]),
						 std::get<0>(date[j][k])
					}
				);
			}
			

		}
	}
}

void BasisSet::createTestBasisSet(const std::vector<Nucleon>& nucleons, double b, double s)
{
	_nucleons = nucleons;
	for (int i = 0; i < nucleons.size(); i++)
	{
		auto date = this->loadAOs(_nucleons[i].charge);
		double tab[] = { 7.402940000 ,  1.576200000 ,  0.3736840000 };

		double c = sqrt(1 / (2 - 2 * exp(-s)));
		date.push_back({
			{s * tab[0], c , {1/sqrt(2*tab[0]),0,0}},
			{s * tab[0], -c, {-1 / sqrt(2 * tab[0]),0,0}}
			});
		date.push_back({
			{s * tab[1], c , {1 / sqrt(2 * tab[1]),0,0}},
			{s * tab[1], -c, {-1 / sqrt(2 * tab[1]),0,0}}
			});
		date.push_back({
			{s * tab[2], c , {1 / sqrt(2 * tab[2]),0,0}},
			{s * tab[2], -c, {-1 / sqrt(2 * tab[2]),0,0}}
			});
		date.push_back({
			{s * tab[0], c , {0,1 / sqrt(2 * tab[0]),0}},
			{s * tab[0], -c, {0,-1 / sqrt(2 * tab[0]),0}}
			});
		date.push_back({
			{s * tab[1], c , {0,1 / sqrt(2 * tab[1]),0}},
			{s * tab[1], -c, {0,-1 / sqrt(2 * tab[1]),0}}
			});
		date.push_back({
			{s * tab[2], c , {0,1 / sqrt(2 * tab[2]),0}},
			{s * tab[2], -c, {0,-1 / sqrt(2 * tab[2]),0}}
			});
		date.push_back({
			{s * tab[0], c , {0,0,1 / sqrt(2 * tab[0])}},
			{s * tab[0], -c, {0,0,-1 / sqrt(2 * tab[0])}}
			});
		date.push_back({
			{s * tab[1], c , {0,0,1 / sqrt(2 * tab[1])}},
			{s * tab[1], -c, {0,0,-1 / sqrt(2 * tab[1])}}
			});
		date.push_back({
			{s * tab[2], c , {0,0,1 / sqrt(2 * tab[2])}},
			{s * tab[2], -c, {0,0,-1 / sqrt(2 * tab[2])}}
			});
		for (int j = 0; j < date.size(); j++)
		{
			_basisFucntion.push_back({});
			for (int k = 0; k < date[j].size(); k++)
			{
				(_basisFucntion.end() - 1)->addPrimitive(
					std::get<1>(date[j][k]), {
						 nucleons[i].p + std::get<2>(date[j][k]),
						 std::get<0>(date[j][k])
					}
				);
			}		
		}

		this->saveAOs(date, _nucleons[i].charge);
	}

}

void BasisSet::calculateOneElectronHamiltonians()
{
	_H_Matrix.resize(_basisFucntion.size());
	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_H_Matrix[i].resize(i+1);
		std::vector<std::future<void>> futures;
		for (int j = 0; j <= i; j++)
		{
			auto wsk = &_H_Matrix[i][j];
			auto bs = &_basisFucntion;
			auto n = &_nucleons;
			futures.push_back(std::async(std::launch::async, [wsk, i, j, bs, n]() {
				*wsk = ((*bs)[i]).calculateHIntegral((*bs)[j], *n); }));
		}
		for (auto& e : futures)
			e.wait();
	}

}

void BasisSet::calculateOverlapIntegrals()
{
	_S_Matrix.resize(_basisFucntion.size());

	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_S_Matrix[i].resize(i + 1);
		std::vector<std::future<void>> futures;
		for (int j = 0; j <= i; j++)
		{
			auto wsk = &_S_Matrix[i][j];
			auto bs = &_basisFucntion;
			futures.push_back(std::async(std::launch::async, [wsk, i, j, bs]() {
				*wsk = (*bs)[i].calculateOverlapIntegral((*bs)[j]); }));
		}
		for (auto& e : futures)
			e.wait();
	}
}

void BasisSet::calculateTwoElectronIntegrals()
{
	_D_Matrix.resize(_basisFucntion.size());
	for (int i = 0; i < _basisFucntion.size(); i++)
	{
		_D_Matrix[i].resize(_basisFucntion.size());
		for (int j = 0; j < _basisFucntion.size(); j++)
		{
			_D_Matrix[i][j].resize(_basisFucntion.size());
			for (int k = 0; k < _basisFucntion.size(); k++)
			{
				_D_Matrix[i][j][k].resize(_basisFucntion.size());
				std::vector<std::future<void>> futures;
				for (int l = 0; l < _basisFucntion.size(); l++)
				{
					auto wsk = &_D_Matrix[i][j][k][l];
					auto bs = &_basisFucntion;
					futures.push_back(std::async(std::launch::async, [wsk, i, j, k, l, bs]() {
						*wsk = ContractedGTO::calulateTwoElectronIntegral((*bs)[i], (*bs)[j], (*bs)[k], (*bs)[l]); }));
				}
				for (auto& e : futures)
					e.wait();
			}
		}

	}
}

void BasisSet::calculateIntegrals()
{
	this->calculateOverlapIntegrals();
	this->calculateOneElectronHamiltonians();
	this->calculateTwoElectronIntegrals();
}

std::vector<std::vector<std::tuple<double, double, Position>>> BasisSet::loadAOs(int n)
{
	std::vector<std::vector<std::tuple<double, double, Position>>> res;
	std::fstream file("data/AO" + std::to_string(n) +".txt", std::ios::in);
	if (!file.is_open())
		throw std::exception("BasisSet::loadAOs couldn't open file");

	int na;
	file >> na;
	res.resize(na);
	for (int i = 0; i < res.size(); i++)
	{
		file >> na;
		res[i].resize(na);
		for (int j = 0; j < res[i].size(); j++)
		{
			double a;
			double b;
			Position p;

			file >> a >> b >> p.x >> p.y >> p.z;
			p.reloadSpherical();
			res[i][j] = std::make_tuple(a, b, p);
		}
	}
	file.close();
	return res;
}

void BasisSet::saveAOs(const std::vector<std::vector<std::tuple<double, double, Position>>>& p, int n)
{
	std::fstream file("data/AO" + std::to_string(n) + ".txt", std::ios::trunc | std::ios::out);
	if (!file.is_open())
		throw std::exception("BasisSet::loadAOs couldn't open file");

	file << (int)p.size();
	file << "\n";
	for (int i = 0; i < p.size(); i++)
	{
		file << p[i].size() << "\n";
		for (int j = 0; j < p[i].size(); j++)
		{
			file << std::get<0>(p[i][j]) << "\t";
			file << std::get<1>(p[i][j]) << "\t";
			file << std::get<2>(p[i][j]).x << "\t";
			file << std::get<2>(p[i][j]).y << "\t";
			file << std::get<2>(p[i][j]).z << "\n";
		}
	}
	file.close();
}

double BasisSet::getSize()
{
	return _basisFucntion.size();
}

double BasisSet::getH(int i, int j)
{
	if (i < j)
		return _H_Matrix[j][i];
	return _H_Matrix[i][j];
}

double BasisSet::getS(int i, int j)
{
	if (i < j)
		return _S_Matrix[j][i];
	return _S_Matrix[i][j];
}

double BasisSet::getDI(int i, int j, int k, int l)
{
	return _D_Matrix[i][j][k][l];
}

double BasisSet::getFunctionValue(int n, const Position& p)
{
	return _basisFucntion[n].F_X(p).real();
}
