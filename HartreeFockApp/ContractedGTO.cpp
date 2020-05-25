#include "ContractedGTO.h"

void ContractedGTO::addPrimitive(double d, const GaussianOrbital& p)
{
	_primitives.push_back({ d,p });
}

double ContractedGTO::calculateOverlapIntegral(const ContractedGTO& right) const
{
	double res = 0.0;
	for (int i = 0; i < this->_primitives.size(); i++)
	{
		for (int j = 0; j < right._primitives.size(); j++)
		{
			res += this->_primitives[i].first * right._primitives[j].first * this->_primitives[i].second.calculateOverlapIntegral(right._primitives[j].second);
		}
	}
	return res;
}

double ContractedGTO::calculateHIntegral(const ContractedGTO& right, const std::vector<Nucleon>& nucleons) const
{
	double res = 0.0;
	for (int i = 0; i < this->_primitives.size(); i++)
	{
		for (int j = 0; j < right._primitives.size(); j++)
		{
			res += this->_primitives[i].first * right._primitives[j].first * this->_primitives[i].second.calculateHIntegral(right._primitives[j].second, nucleons);
		}
	}
	return res;
}

double ContractedGTO::calculateHIntegral(const ContractedGTO& right, const std::vector<Nucleon>& nucleons, double S) const		//TODO
{
	return 0.0;
}

double ContractedGTO::calulateTwoElectronIntegral(const ContractedGTO& a, const ContractedGTO& b, const ContractedGTO& c, const ContractedGTO& d)
{
	double res = 0.0;
	for (int i = 0; i < a._primitives.size(); i++)
	{
		for (int j = 0; j < b._primitives.size(); j++)
		{
			for (int k = 0; k < c._primitives.size(); k++)
			{
				for (int l = 0; l < d._primitives.size(); l++)
				{
					res += a._primitives[i].first * b._primitives[j].first * c._primitives[k].first * d._primitives[l].first * 
						GaussianOrbital::calulateTwoElectronIntegral(a._primitives[i].second, b._primitives[j].second, 
							c._primitives[k].second, d._primitives[l].second);
				}
			}
		}
	}

	return res;
}

double ContractedGTO::calulateTwoElectronIntegral(const ContractedGTO& a, const ContractedGTO& b, const ContractedGTO& c, const ContractedGTO& d, double Spq, double Srs)//TODO
{
	return 0.0;
}

std::complex<double> ContractedGTO::F_X(const Position&p)
{
	std::complex<double> res = { 0,0 };
	for (int i = 0; i < _primitives.size(); i++)
	{
		res += std::complex<double>{_primitives[i].first}*_primitives[i].second.F_X(p);
	}
	return res;
}
