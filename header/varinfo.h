#ifndef __VAR_INFO_H__
#define __VAR_INFO_H__

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <pthread.h>
#include <stdexcept>
#include <iomanip>
#include <tr1/unordered_map>
#include "stdlib.h"
#include <algorithm>
#include <numeric>
#include <ctime>
#include <time.h>

#include "./alglib.h"
#include "./bam/bam.h"
#include "./bio_files.h"
#include "./util.h"
#include "./suffix_array.h"


#define EXTRA_LEN	500

using namespace std;

class CINFO;

struct QCINFO
{
	float fVaf;
	float fStrandBias;
	vector<uint16_t> vnReadLen;
	vector<uint8_t> vnMapQ;
	vector<float> vfBaseQ;
	uint16_t nFrontA;
	uint16_t nFrontC;
	uint16_t nFrontG;
	uint16_t nFrontT;
	uint16_t nReverseA;
	uint16_t nReverseC;
	uint16_t nReverseG;
	uint16_t nReverseT;
	vector<string> vsForwardIns;
	vector<uint16_t> vnForwardDel;
	vector<string> vsReverseIns;
	vector<uint16_t> vnReverseDel;
	string sContingencyTable;
	QCINFO(){
		sContingencyTable = "0,0,0,0";
		nFrontA = nFrontC = nFrontG = nFrontT = 0;
		nReverseA = nReverseC = nReverseG = nReverseT = 0;
	}
};

struct VARIANT
{
	vector<bool> vbIsPass;

	//raw data
	vector<string> vsHeader;
	
	vector<bool> vbIsVcf;
	vector<string> vsChr;
	vector<int>	vnPos;
	vector<string> vsRef;
	vector<string> vsAlt;
	vector<double> vdVaf;
	
	vector<string> vsRaw;

	//AS IS: BAM QC Information
	// vector<float> vfVaf;					// variant allele freq
	vector<string> vsContingencyTable;
	vector<float> vfStrandBias;				// strand bias 
	vector<vector<uint16_t> > v2nReadLen;	// [loci][read lenths]
	vector<vector<uint8_t> > v2nMapQ;		// [loci][mapping qualities]
	vector<vector<float> > v2fBaseQ;		// [loci][base qualities]
};

struct DIST_THREAD
{
	CINFO *INFO;
	int nId;
};



class CINFO
{
public:
	bool m_bIsDebug;
	string m_sBamFile;
	bool m_bIsMod;
	string m_sInputFile;
	string m_sOutputFile;
	int m_nCntThread;
	double m_dVaf;


	//for util
	clock_t m_clockStart;
	clock_t m_clockEnd;
	struct timespec m_tspecStart;
	struct timespec m_tspecEnd;
public:
	CINFO(string, string, string, bool, int, double, bool);
	bool ReadInput();
	bool ModVariant();
	bool CalcInfo();
	bool Report();
	
	int PrintCommonInfo();

private:
	// input
	VARIANT m_Input;
	
	// calc info
	bool DelIns(int, int);
	string GetRead(bam1_t*, int);

	bool GetAnalysis(QCINFO &,string,string);
	bool GetAnalysisDelIns(QCINFO &, string, string);
	bool GetNb(bam1_t *, int, QCINFO &);
	bool GetDelIns(bam1_t *, int, QCINFO &);
	void* AlleleDist(int);
	static void *AlleleDist_helper(void *object)
	{
		DIST_THREAD *my = (DIST_THREAD*)object;
		my->INFO->AlleleDist(my->nId);
	}
	bool GetVarInfo(int, string, int, int, string, string, bam_index_t *, bamFile);
	int ConvertChrToTid(string, bam_header_t*);
/*
	bool GetReadAlign(int &, int &, string &, int, string, string, bam1_t*);
	bool CheckRepeat(string, string, string, int, int &, int &);
	bool AlleleCount(bam_index_t *, bamFile &, int, int, string, string, int, int, int &, vector<int> &, vector<int> &, vector<int> &, vector<int> &, vector<string> &, vector<string> &);
	//double CalcCoef(string sChr, vector<string> &, vector<int> &, vector<int> &);
	bool FilterRead(string, vector<int> &, vector<int> &, vector<int> &, vector<int> &, vector<string> &, vector<string> &, int &, int &, int &, int &);
*/
	// report 
	bool ReportBamInfo();
//	bool ReportAdiscan();
//	bool ReportVarscan();


	// library
	bool Parsing(vector<string> &, string, string);
	bool CalcPearsonCorr(int*, int*, double&, double&);

	//util
	void SetStartTime()
	{
		//m_clockStart = clock();
		clock_gettime(CLOCK_REALTIME, &m_tspecStart);
	};
	int ViewStatus(int nCurr, int nTotal, string sId, bool bIsNoDel=false);
	string to_string(int);
	string to_string(unsigned long);
	string to_string(uint16_t);
	string to_string(double);
	string Table2String(int,int,int, int);
};





#endif

