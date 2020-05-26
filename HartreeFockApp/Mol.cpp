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
	_orbitalCoeficients[0].resize(_basisSet.getSize(), c);
	for (int i = 0; i < c; i++)
	{
		for (int j = 0; j < _basisSet.getSize(); j++)
		{
			_orbitalCoeficients[0](j, i) = 0.0;
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
		for (int i = 4; i > 0; i--)
		{
			_orbitalCoeficients[i] = _orbitalCoeficients[i - 1];
		}
		this->recalculateP();
		this->recalculateF();
		EigenSolver::solve(_FMatrix, _SMatrix, _orbitalCoeficients[0]);


		//this->PullyMixing(std::min(n,5));
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
		if (i > 300)
			break;

	} while (abs(oldE - newE) > delta || i < 10);
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
		res[i] = _orbitalCoeficients[0](i,m);
	}
	return res;
}

double Mol::countMolecularFunction(int m, const Position & p)
{
	double res = 0.0;
	for (int i = 0; i < _basisSet.getSize(); i++)
	{
		res += _basisSet.getFunctionValue(i, p) * _orbitalCoeficients[0](i, m);
	}
	return res;
}

double Mol::countMolecularFunction(int m, double x)
{
	return this->countMolecularFunction(m, Position{ x,0.0,0.0 });
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
				s += _orbitalCoeficients[0](i, k) * _orbitalCoeficients[0](j, k);
			}
			_PMatrix(i, j) =  2 * s;
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

void Mol::PullyMixing(int n)
{
	if (n == 0)
		return;
	std::vector<arma::mat> B(_orbitalCoeficients[0].n_cols);
	std::vector<arma::vec> X(_orbitalCoeficients[0].n_cols);
	arma::vec L(n);
	for (int i = 0; i < n-1; i++)
	{
		L(i) = 0;
	}
	L(n - 1) = -1;
		
	for (size_t i = 0; i <_orbitalCoeficients[0].n_cols; i++)
	{
		B[i] = arma::mat(n, n);
		X[i] = arma::vec(n);
	}

	for (int i = 0; i < n-1; i++)
	{
		for (int j = 0; j < n-1; j++)
		{
			auto vec = this->getCoeficientIntegrals(_orbitalCoeficients[i] - _orbitalCoeficients[i + 1], _orbitalCoeficients[j] - _orbitalCoeficients[j + 1]);
			for (size_t k = 0; k < _orbitalCoeficients[0].n_cols; k++)
			{
				B[k](i, j) = vec[k];
			}
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < B.size(); j++)
		{
			B[j](i, n-1) = -1;
			B[j](n-1, i) = -1;
			B[j](n-1, n-1) = 0;
		}		
	}
	for (int i = 0; i < B.size(); i++)
	{
		X[i] = arma::solve(B[i], L);
	}

	for (int i = 0; i < _orbitalCoeficients[0].n_cols; i++)
	{
		for (int j = 0; j < _orbitalCoeficients[0].n_rows; j++)
		{

			_orbitalCoeficients[0](j, i) /=2;
			for (int k = 0; k < n-1; k++)
			{
				_orbitalCoeficients[0](j, i) += 0.5 * X[i](k) * _orbitalCoeficients[k + 1](j, i);
			}
		}
	}

	for (int i = 0; i < _orbitalCoeficients[0].n_cols; i++)
	{
		double res = 0.0;
		for (int j = 0; j < _orbitalCoeficients[0].n_rows; j++)
		{
			for (int k = 0; k < _orbitalCoeficients[0].n_rows; k++)
			{
				res += _orbitalCoeficients[0](j, i) * _orbitalCoeficients[0](k, i) * _SMatrix(k, j);
			}
		}
		double l = sqrt(res);
		for (int j = 0; j < _orbitalCoeficients[0].n_rows; j++)
			_orbitalCoeficients[0](j, i) /= l;
	}
}


