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
	vector<int> nReadLen;
	vector<uint8_t> nMapQ; // mapping quality
	vector<uint8_t> nBaseQ; // base quality
	uint16_t nReverse;
	uint16_t nDepthA;
	uint16_t nDepthC;
	uint16_t nDepthG;
	uint16_t nDepthT;
	uint16_t nDepthIndel;
	QCINFO(){
		nReverse = nDepthA = nDepthC = nDepthG = nDepthT = nDepthIndel = 0;
	}
};

struct VARIANT
{
	//raw data
	vector<string> vsHeader;
	
	vector<bool> vbIsVcf;
	vector<string> vsChr;
	vector<int>	vnPos;
	vector<string> vsRef;
	vector<string> vsAlt;
	
	vector<string> vsRaw;



	//TODO:information
	QCINFO QcInfo;

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

	//RDscan m_rdscan;			//rdscan (output)

	//for util
	clock_t m_clockStart;
	clock_t m_clockEnd;
	struct timespec m_tspecStart;
	struct timespec m_tspecEnd;
public:
	CINFO(string, string, string, bool, int, bool);
	bool ReadInput();
	bool ModVariant();
	bool CalcInfo();
	//bool Report();
	
	int PrintCommonInfo();

private:
	// input
	VARIANT m_Input;
	
	// calc info
	bool GetStrandBias(QCINFO &);
	bool GetVaf(QCINFO &);
	bool GetNb(bam1_t *, int, QCINFO &);
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

	// report 
	bool ReportVcf();
//	bool ReportAdiscan();
//	bool ReportVarscan();

*/
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

};





#endif

