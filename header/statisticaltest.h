#ifndef STATISTICAL_TEST_H_
#define STATISTICAL_TEST_H_

#include <iostream>
#include <cmath>

namespace STATTEST
{
	//================== << Fisher's Exact Test >> ==================//
	/*
	**    n11  n12  | n1_
	**    n21  n22  | n2_
	**   -----------+----
	**    n_1  n_2  | n
	**/
	double GetFisherPvalue(int n11, int n12, int n21, int n22, bool twosided=false);
	double LogHypergeometricProb(int, int, int, int);
	double LogFac(int n);
	//===============================================================//
}

#endif