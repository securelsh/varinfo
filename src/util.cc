#include "../header/varinfo.h"
pthread_mutex_t mutex_util = PTHREAD_MUTEX_INITIALIZER;

int CINFO::PrintCommonInfo()
{
	cout << "####################################" << endl;
	cout << "## VARinfo v0.0.1 (31 Mar 2021)   ##" << endl;
	cout << "## (C) 2021-2022 Syntekabio Inc.  ##" << endl;
	cout << "####################################" << endl;

	cout << endl;
	cout << "[INPUT Files]" << endl;
	cout << "- Input File   : " << m_sInputFile << endl;
	cout << "- Output File  : " << m_sOutputFile << endl;
	cout << "- Bam File     : " << m_sBamFile << endl;
	cout << "[OPTIONs]" << endl;
	cout << "- # of threads : " << m_nCntThread;
	if(m_bIsMod)	cout << "- Consecutive SNP to Indel" << endl;

	return 0;
}

/*
bool CRD::Parsing(vector<string> &vsWord, string sLine, string sSep)
{
	char *pcEnd;
	unsigned long ulMarker;
	vsWord.clear();
	while(1)
	{
		ulMarker = sLine.find_first_of(sSep);
		vsWord.push_back(sLine.substr(0, ulMarker));

		if(ulMarker > sLine.size()) break;
		sLine = sLine.substr(ulMarker+1);
	}
	return true;
}


string CRD::to_string(uint16_t nNum)
{
	stringstream ss;
	ss << nNum;
	return ss.str();
}


string CRD::to_string(int nNum)
{
	stringstream ss;
	ss << nNum;
	return ss.str();
}

string CRD::to_string(unsigned long ulNum)
{
	stringstream ss;
	ss << ulNum;
	return ss.str();
}

string CRD::to_string(double dNum)
{
	stringstream ss;
	ss << dNum;
	return ss.str();
}

//////////////////////////////////////////// Loading view ////////////////////////////////
int CRD::ViewStatus(int nCurr, int nTotal, string sId, bool bIsNoDel)
{
	//m_clockEnd = clock();
	clock_gettime(CLOCK_REALTIME, &m_tspecEnd);
	
	pthread_mutex_lock(&mutex_util);
	printf("\r[PROCESS]: %s (%d/%d) [", sId.c_str(), nCurr, nTotal);

	double dPercent = (double)(nCurr)/(double)nTotal*100;
	int nBarCount = dPercent/2;
	for(int k=0; k<50; k++){
		if(nBarCount > k)       printf("%c", '=');
		else                    printf("%c", ' ');
	}
	int nSec = m_tspecEnd.tv_sec - m_tspecStart.tv_sec;
	int nMin = 0;
	int nHour = 0;
	
	nHour = nSec/3600;
	nSec -= 3600*nHour;
	nMin = nSec/60;
	nSec -= 60*nMin;

	if(nHour > 0)
	{
		if(bIsNoDel)	printf("] %0.2f%%  (%dh %dm %ds)           \n", dPercent, nHour, nMin, nSec);
		else			printf("] %0.2f%%  (%dh %dm %ds)            ", dPercent, nHour, nMin, nSec); //(int)((m_clockEnd-m_clockStart)/CLOCKS_PER_SEC));
	}
	else if(nMin > 0)
	{
		if(bIsNoDel)	printf("] %0.2f%%  (%dm %ds)               \n", dPercent, nMin, nSec);
		else			printf("] %0.2f%%  (%dm %ds)                ", dPercent, nMin, nSec); //(int)((m_clockEnd-m_clockStart)/CLOCKS_PER_SEC));
	}
	else
	{
		if(bIsNoDel)	printf("] %0.2f%%  (%ds)                   \n", dPercent, nSec);
		else			printf("] %0.2f%%  (%ds)                    ", dPercent, nSec); //(int)((m_clockEnd-m_clockStart)/CLOCKS_PER_SEC));
	}
	
	
	
	fflush(stdout);
	pthread_mutex_unlock(&mutex_util);

	return 0;
}
*/
