#ifndef __PROGPARAMS_H__
#define __PROGPARAMS_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <math.h>

#define DEBUG 0
#define DBL_MAX  144115188075855872

// a class to read and manage parameters

class programParams
{

public:
	programParams();
	int initialize(char *filename);
	void announce();

	int matrixsize;
	int stepnumber;
	char *outfile;
	int initialconfig;
	char *inifile;
	double Lc;
	double T0;
	double Tf;
	double kfac;
	double wmoveA;
	int Type;

};


// This is the end of the header guard
#endif
