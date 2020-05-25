#pragma once
#include <armadillo>
#include <chrono>
#include <future>

class EigenSolver
{
public:
	static void solve(const arma::mat& H, const arma::mat & S, arma::mat & coM)
	{
		auto coefs = ortogonalizeSet(S);
		/*for (size_t i = 0; i < coefs.n_rows; i++)
		{
			for (size_t j = 0; j < coefs.n_cols; j++)
			{
				std::cout << std::ceil(coefs(i,j) * 100) / 100.0 << "\t";
			}
			std::cout << std::endl;
		}*/
		//std::cout << std::endl;
		arma::mat _H_ = H;
		for (int i = 0; i < H.n_rows; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				double a = countMatrixElement(i, j, coefs, H);
				_H_(i, j) = a;
				_H_(j, i) = _H_(i, j);
			}

		}
		arma::vec eigval;
		arma::mat eigvec;	
		arma::eig_sym(eigval, eigvec, _H_);

		for (int i = 0; i < coM.n_cols; i++)
		{
			double test = 0, set = 0;
			for (int j = 0; j < coM.n_rows; j++)
			{
				double s = 0.0;
				for (int k = 0; k < eigvec.n_rows; k++)
				{
					s += eigvec(k, i) * coefs(k, j);
				}
				coM(j, i) = s;
			}
			//std::cout<<eigval[0] << std::endl;
		}
		/*for (size_t i = 0; i < eigval.size(); i++)
		{
			std::cout << std::ceil(eigval[i] * 1000) / 1000.0 << "\t";
		}
		std::cout << std::endl;*/
	}


	static arma::mat ortogonalizeSet(const arma::mat& S)
	{
		int n = S.n_rows;
		arma::mat res(S.n_rows, S.n_cols);

		for (int a = 0; a < n; a++)
		{
			for (int i = 0; i < n; i++)
			{
				if (i < a)
				{
					double sum = 0.0;
					for (int k = 0; k < a; k++)
					{
						for (int l = 0; l < a; l++)
						{
							sum -= res(k, i) * res(k, l) * S(l, a);
						}
					}
					res(a, i) = sum;
				}
				else if (i == a)
					res(a, i) = 1.0;
				else
					res(a, i) = 0;
			}

			double l = sqrt(countMatrixElement(a, a, res, S));
			for (int i = 0; i < n; i++)
			{
				res(a, i) /= l;
			}

		}
		return res;
	}

	static double countMatrixElement(int r, int c, const arma::mat& coeficcient, const arma::mat& M)
	{
		double res = 0.0;
		for (int i = 0; i < coeficcient.n_cols; i++)
		{
			for (int j = 0; j < coeficcient.n_cols; j++)
			{
				res += coeficcient(r, i) * coeficcient(c, j) * M(i, j);
			}
		}
		return res;
	}

};