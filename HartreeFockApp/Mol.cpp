#include "Mol.h"

Mol::Mol(const std::vector<Nucleon>& n)
{
	_nucleons = n;
	for (int i = 0; i < n.size(); i++)
	{
		_MOcount += n[i].charge;
	}
	_MOcount = (_MOcount + 1) / 2;
	setMOcount(_MOcount);

}

Mol::Mol()
{
}

Mol::Mol(int charge):Mol({ {Position{0,0,0}, charge } })
{
}

void Mol::addNucleon(const Nucleon& n)
{
	_nucleons.push_back(n);
	_MOcount = 0;
	for (int i = 0; i < _nucleons.size(); i++)
	{
		_MOcount += _nucleons[i].charge;
	}
	_MOcount = (_MOcount + 1) / 2;
	setMOcount(_MOcount);
}

void Mol::moveNucleon(int n, const Position& delta)
{
	_nucleons[n].p = _nucleons[n].p + delta;
	_integralsCalculated = false;
	_setInitiated = false;
}

void Mol::setMOcount(int c)
{
	_MOcount = c;
	_orbitalCoeficients.resize(_basisSet.getSize(), c);
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < _basisSet.getSize(); j++)
		{
			_orbitalCoeficients(j, i) = 0.0;
		}
	}
}

void Mol::initBasisSet()
{
	_basisSet.createBasisSet(_nucleons);
	this->setMOcount(_MOcount);
	_setInitiated = true;
	_FMatrix.resize(_basisSet.getSize(), _basisSet.getSize());
	_SMatrix.resize(_basisSet.getSize(), _basisSet.getSize());
	_PMatrix.resize(_basisSet.getSize(), _basisSet.getSize());
}

void Mol::calculateIntegrals()
{
	if (_setInitiated)
	{
		_basisSet.calculateIntegrals();
		_integralsCalculated = true;
		for (int r = 0; r < _basisSet.getSize(); r++)
		{
			for (int s = 0; s <= r; s++)
			{
				_SMatrix(r, s) = _basisSet.getS(r, s);
				_SMatrix(s, r) = _basisSet.getS(r, s);
			}
		}
	}
		
}

void Mol::HFProcedure(int n)
{
	if (_integralsCalculated)
	{
		_oldPMatrix = _PMatrix;
		EigenSolver::solve(_FMatrix, _SMatrix, _orbitalCoeficients);
		this->recalculateP();
		this->recalculateF();
		this->recalculateEnergy();

	}
}

void Mol::HF_TO_Divergance(double delta)
{
	double oldE;
	double newE = 1000000;
	int i = 0;
	do {
		HFProcedure(i);
		oldE = newE;
		newE = this->getElectronicEnergy();
		i++;
		if (i > 150)
		{
			std::cout << "Energy doesn't converge" << std::endl;
			break;
		}
			

	} while (abs(oldE-newE) > delta);

}

void Mol::HF_TO_ElapsedTime(double t)
{
}

double Mol::getElectronicEnergy()
{
	return _electronicEnergy;
}

double Mol::getMoleculeEnergy()
{
	auto nE = 0.0;
	for (int i = 0; i < _nucleons.size(); i++)
	{
		for (int j = 0; j < i; j++)
		{
			nE += _nucleons[i].charge * _nucleons[j].charge / (_nucleons[i].p - _nucleons[j].p).r;
		}
	}

	return nE + _electronicEnergy;
}

std::vector<double> Mol::getMolecularCoeficents(int m)
{
	std::vector<double> res;
	res.resize(_basisSet.getSize());

	for (int i = 0; i < _basisSet.getSize(); i++)
	{
		res[i] = _orbitalCoeficients(i,m);
	}
	return res;
}

double Mol::countMolecularFunction(int m, const Position & p)
{
	double res = 0.0;
	for (int i = 0; i < _basisSet.getSize(); i++)
	{
		res += _basisSet.getFunctionValue(i, p) * _orbitalCoeficients(i, m);
	}
	return res;
}

double Mol::countMolecularFunction(int m, double x)
{
	return this->countMolecularFunction(m, Position{ x,0.0,0.0 });
}

void Mol::saveMolecularFunctionPlane(const std::string& filename, double yMax, double zMax, double delta)
{
	for (int i = 0; i < _orbitalCoeficients.n_cols; i++)
	{
		std::fstream file("data/" + filename +"_"+ std::to_string(i) + ".txt", std::ios::trunc | std::ios::out);
		if (!file.is_open())
			throw std::exception("Mol::saveMolecularFunctionPlane couldn't create file");

		double z = -zMax;
		while (z <= zMax)
		{
			double y = -yMax;
			while (y <= yMax)
			{
				auto v = this->countMolecularFunction(i, { 0.0,y,z });
				file << v <<"\n";

				y += delta;
			}
			z += delta;
		}

		file.close();
	}
}

void Mol::recalculateP()
{
	for (int i = 0; i < _basisSet.getSize(); i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double s = 0.0;
			for (int k = 0; k < _MOcount; k++)
			{
				s += _orbitalCoeficients(i, k) * _orbitalCoeficients(j, k);
			}
			_PMatrix(i, j) =  2 * s;
			_PMatrix(i, j) = _PMatrix(i, j) * 0.5 + _oldPMatrix(i, j) * 0.5;
			_PMatrix(j, i) = _PMatrix(i, j);

		}
	}

}

void Mol::recalculateF()
{
	for (int r = 0; r < _basisSet.getSize(); r++)
	{
		for (int s = 0; s <= r; s++)
		{
			double sum = _basisSet.getH(r, s);
			for (int p = 0; p < _basisSet.getSize(); p++)
			{
				for (int q = 0; q < _basisSet.getSize(); q++)
				{
					sum += _PMatrix(q, p) * (_basisSet.getDI(r, p, s, q) - _basisSet.getDI(r, p, q, s) / 2.0);
				}
			}
			_FMatrix(r, s) = sum;
			_FMatrix(s, r) = sum;
		}
	}
}

void Mol::recalculateEnergy()
{
	_electronicEnergy = 0.0;
	for (int r = 0; r < _basisSet.getSize(); r++)
	{
		for (int s = 0; s < _basisSet.getSize(); s++)
		{
			_electronicEnergy += _PMatrix(s, r) * (_basisSet.getH(s, r) + _FMatrix(r, s));
		}
	}
	_electronicEnergy /= 2.0;
}

arma::vec Mol::getCoeficientIntegrals(const arma::mat& left, const arma::mat& right)
{
	auto res = arma::vec(left.n_cols);

	for (int i = 0; i < left.n_cols; i++)
	{
		res(i) = 0.0;
		for (int j = 0; j < left.n_rows; j++)
		{
			for (int k = 0; k< left.n_rows; k++)
			{
				res(i) += left(j, i) * right(k, i) * _SMatrix(j, k);
			}
		}
	}
	return res;
}


bool Mol::hardDivergance(double d)
{
	for (int i = 0; i < _orbitalCoeficients.n_rows; i++)
	{
		for (int j = 0; j < _orbitalCoeficients.n_cols; j++)
		{
			if (abs(_orbitalCoeficients(i, j) - _orbitalCoeficients(i, j)) > d)
				return true;
		}
	}
	return false;
}


