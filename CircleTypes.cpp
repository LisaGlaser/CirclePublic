#include "CircleTypes.h"

std::complex<double> FixedRealCircle::JDCoefficients (int i, int j)
{
	std::complex<double> c;
	int sign;
	if ((i != j) && (i != -j)){
		c = 0;
	}
	else if (i>0){
		c = 1 / sqrt(2);
	}
	else if (i<0){
		if (j > 0){
			sign = -1;
		}
		else{
			sign = 1;
		}
		std::complex<double> I(0.0, 1.0);
		c = (std::complex<double>)(sign) * I / sqrt(2);
	}
	else if (i == 0){
		c = 1;
	}
	return c;
}

void FixedRealCircle::makeJDiagonalizer()
{
	JDiagonalizer = Eigen::MatrixXcd(size, size);
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){
			int possize = int((size - 1) / 2);
			int m = possize - i; //pluses come first, zero in the middle
			int n = possize - j;
			JDiagonalizer(i, j) = JDCoefficients(m, n);
		}
	}
}

Eigen::MatrixXcd FixedRealCircle::getJDiagonalizer(){
	return JDiagonalizer;
}

void FixedRealCircle::setSize(int newsize){
	size = 2*newsize+1; // ensure odd sizes
	makeJDiagonalizer();
}

int FixedRealCircle::getSize(){
	return size;
}
