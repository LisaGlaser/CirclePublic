#include "Dirac.h"


Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es;
const std::complex<double> imagUnit(0.0, 1.0);

Dirac::Dirac(programParams iniV)
{
	type=iniV.Type;
	size=iniV.matrixsize;
	FixedRealCircle Circle;
	Circle.setSize(size);
	size = Circle.getSize();
	truesize = size;
	JDiagonalizer = Circle.getJDiagonalizer();

  	// the dimension of D is twice the argument supplied here, so 2*size
	setU();
	initial(iniV.initialconfig,iniV.inifile);
	makeD();
	es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	LcI=iniV.Lc; // I eventually want to remove Lc from the code completely
	Lc=LcI;
	LcFactors();
	calc_RVd();
	calc_S1con();
	moved=0;

}

Dirac::~Dirac()
{
	char finale[250];
	std::cout<<"Deconstructing Dirac"<<std::endl;
	strcpy(finale,"deconstructor_finalmatrix.txt");
	//printAll(finale);

	// This desctructor should also clean up my entire memory
	/// and I need to catch the sigkill signal to make it work
	//delete D;

}

void Dirac::dynamicLc()
{
	// no need to initialize, any square is larger than a negative number
	// unless we get imaginary numbers here, but then we are screwed anyway
	for(int i=0;i<eivals.size();i++)
	{
		if(eivals[i]*eivals[i]>Lc) Lc=eivals[i]*eivals[i];
	}
 //number from Mathematica
}

void Dirac::LcFactors()
{
	if(LcI<0) dynamicLc();
	E10=2.*log(Lc+1)/(Lc+.001);
	// numbers fixed on 06/06/18 hopefully now correct
	fLd=-41.928230072331301668/sqrt(E10); // number copied from Mathematica
	volumecoeff=4.67745*sqrt(E10);
}




double Dirac::F10(double l)
{
  return exp(-l*E10)/(1. + l*E10);
}

double Dirac::F11(double l)
{
	double q1=0.70710678118654752440;
	return exp(-1.-l*E10)/(1. + l*E10) - q1*exp(-1. - 2.*l*E10)/(0.5 + l*E10);
}


void Dirac::calc_RVd()
{
	/// calculates R and V both because they go together so well
	double tunedM = 0.39612512715394627; // m is tuned so that ds = 2.0 for the sphere at size 30
	double S=0,rpart=0,volpart=0,kt=0,ktl2=0,ex=0;
	double ep=tunedM*log(Lc+1)*E10; // where 0.15 is the abel special factor

	// We got the normal eigenvalues and can just square them
	// and then throw away eigenvalues larger than Lc

	if(LcI<0) LcFactors(); /// actually not so good, now it does not recompute

	for(int x=0;x<eivals.size();x++)
	{
		ex=eivals[x];

		rpart+=F11(ex*ex);
	    volpart+=F10(ex*ex);

		kt+=exp(-1.*ex*ex*ep);
		ktl2+=ex*ex*exp(-1.*ex*ex*ep);

	}

	ds=2*ep*ktl2/kt; // get the spectral dimension by evaluating t p(t,D^2) / p(t,D) at t=\epsilon
	R=fLd*rpart;
	Vol=volumecoeff*volpart;
}

void Dirac::calc_S1con()
{
		std::complex <double> S1c=0;
		Eigen::MatrixXcd tempM1(truesize,truesize);
		tempM1=D*U-U*D;
		tempM1=U.adjoint()*tempM1;
		tempM1=tempM1-Eigen::MatrixXcd::Identity(truesize, truesize);
		Eigen::Map<Eigen::RowVectorXcd> v1(tempM1.data(), tempM1.size());
		S1c=v1*v1.adjoint();

		con=S1c.real();
	}


void Dirac::initial(int ini,char *inifile)
{

	Eigen::MatrixXcd mtemp(size,size);

	if(ini == 1) // a real circle spectrum in H0
	{
		Eigen::MatrixXcd stdD = Eigen::MatrixXcd::Zero(size, size);
		for (int i = 0; i < size; i++){
			for (int j = 0; j < size; j++){
				int possize = int((size - 1) / 2);
				int n = possize - i; //pluses come first, zero in the middle
				int m = possize - j;
				stdD(i, i) = n;
			}
		}

		Eigen::MatrixXcd mtemp(size,size);
		mtemp = JDiagonalizer * stdD * JDiagonalizer.adjoint(); // purely imaginary, hence anticommuting with J
		H0 = mtemp.imag();
	}
	else if(ini == 2) // a random spectrum in H0
	{
		Eigen::MatrixXcd stdD = Eigen::MatrixXd::Random(size,size);


		Eigen::MatrixXcd mtemp(size,size);
		mtemp = JDiagonalizer * stdD * JDiagonalizer.adjoint(); // purely imaginary, hence anticommuting with J
		H0 = mtemp.imag();
	}
	else if(ini == 3) // loading from an initial condition
	{
		H0.resize(size,size);

		std::vector<std::complex<double>> vecH0;
		std::vector<double> vecH0r;
		HighFive::File file(inifile, HighFive::File::ReadOnly);
		std::string DATASET_NAME("WTF?");
		if(type==0)
		{
			DATASET_NAME="H0";
		}
		else if(type==2)
		{
			DATASET_NAME="H0, L0";
		}

		// we get the dataset
		HighFive::DataSet dataset = file.getDataSet(DATASET_NAME);

		// we convert the hdf5 dataset to a single dimension vector
		dataset.read(vecH0);

		for(int i=0; i<vecH0.size(); i++)
		{
 		vecH0r.push_back(vecH0[i].real());
		}

		read_H0(vecH0r);
	}
	else
	{
		std::cout<< " Sorry we don't support that many initial configurations yet"<<std::endl;
	}


}

std::vector<std::complex<double>> Dirac::get_D()
{

	std::vector<std::complex<double>> Dvec;

	for(int x=0;x<truesize;x++)
	{
		for(int y=0;y<truesize;y++)
		{
			Dvec.push_back(D(x,y));
		}
	}
	return Dvec;
}

std::vector<double> Dirac::get_final()
{

	std::vector<double> Dvec;

		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				Dvec.push_back(H0(x,y));
			}
		}

	return Dvec;
}

void Dirac::read_H0(std::vector<double> H0v)
{

		for(int x=0;x<size;x++)
		{
			for(int y=0;y<size;y++)
			{
				//Dvec.push_back(H0(x,y));
				H0(x,y)=H0v[x*size+y];
			}
		}

}

std::vector <double> Dirac::get_eigenvalues()
{
	std::vector <double> evs;

	for(int x=0;x<eivals.size();x++)
	{
	 	evs.push_back(eivals[x]);

	 }

	return evs;
}

double Dirac::getR()
{
	if(moved==1)
	{
		es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	calc_S1con();
	calc_RVd();
	}
	moved=0;
	return R;
}

double Dirac::getVol()
{
	if(moved==1)
	{
		es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	calc_S1con();
	calc_RVd();
	}
	moved=0;
	return Vol;
}


double Dirac::getds()
{
	if(moved==1)
	{
		es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	calc_S1con();
	calc_RVd();
	}
	moved=0;
	return ds;
}


void Dirac::setU()
{
	//std::cout<<truesize<<std::endl;
	U = Eigen::MatrixXcd::Zero(truesize,truesize);
	U.topRightCorner(truesize-1,truesize-1) = Eigen::MatrixXcd::Identity(truesize-1,truesize-1);

}


double Dirac::getcon()
{
	if(moved==1)
	{
	es.compute(D,false); // computes the eigenvalues, false tells it not to compute the eigenvectors
	eivals=es.eigenvalues();
	calc_S1con();
	calc_RVd();
	}
	moved=0;
	return con;
}


void Dirac::moveA(double pA)
{
	int hermitian;
	Eigen::MatrixXd temp;

	hermitian=-1;
	moveA_raw(H0,hermitian,pA,size);

	if(!H0.isApprox(hermitian*H0.adjoint())) std::cout<<"H0 warning"<<std::endl;
	temp=H0+hermitian*H0.adjoint();
	H0=0.5*temp;

	makeD();
	moved=1;
	if(!D.isApprox(D.adjoint()))
	{
		std::cout<<"That Dirac is not selfadjoint after moveM, something is wrong here!"<<std::endl;
		std::cout<<D-D.adjoint()<<std::endl;
	}

}



int seed = (int)time(0);    // random seed
StochasticLib1 sto(seed);           // make instance of random library



int moveA_raw(Eigen::MatrixXd &M, int her, double p,int s)
{
	Eigen::MatrixXd diff(s,s);
	diff=Eigen::MatrixXd::Random(s,s);

	M+=sto.Normal(p,p)*(diff+ her*diff.adjoint()); // here comes the anti hermitian


	return 0;

}

void Dirac::makeD()
{
	// H0 is antisymmetric.
	// we have <f_i, D f_j> = - i (H0)_{ij}
	Eigen::MatrixXcd H = imagUnit * H0.cast<std::complex<double>>();
	D = JDiagonalizer.adjoint() * H * JDiagonalizer;
}
