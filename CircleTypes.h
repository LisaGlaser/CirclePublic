#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include <iostream>

// creates the parametrisations of the dirac operator needed for the code.

class FixedRealCircle
{
protected:
	int size;
	Eigen::MatrixXcd JDiagonalizer; // this is the matrix U s.t. J f_n = f_n with e_n = U f_n, i.e. U_{ij} = <f_i, e_j>
	void makeJDiagonalizer();
	std::complex<double> JDCoefficients (int i, int j);
public:
	int getSize();
	void setSize(int newsize);
	Eigen::MatrixXcd getJDiagonalizer();
};
