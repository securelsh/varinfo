#ifndef STATISTICAL_TEST_H_
#define STATISTICAL_TEST_H_

#include <iostream>
#include <cmath>

/*
 - Changelog:
   -- version 0.1.0 Seok Chol Hong shulkhorn@gmail.com
      --- Fisher's exact test with logarithmic scale for large values
*/
namespace STATTEST
{
	double LogFac(int);
	void InitLogFacs(double*, int);

	//================== << Fisher's Exact Test >> ==================//
	/*
	**    n11  n12  | n1_
	**    n21  n22  | n2_
	**   -----------+----
	**    n_1  n_2  | n
	**
	** return two-tailed P-value
	**/
	double GetFisherPvalue(int n11, int n12, int n21, int n22);
	double FastGetFisherPvalue(int n11, int n12, int n21, int n22);
	double LogHypergeometricProb(int, int, int, int);
	double FastLogHypergeometricProb(double*, int, int, int, int);
	//===============================================================//
}

#endif