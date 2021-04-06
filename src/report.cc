#include "../header/varinfo.h"


using namespace std;


bool CINFO::Report()
{
	// output mod adiscan
	ofstream fout;
	fout.open(m_sOutputFile.c_str());
	for(int i=0; i<(int)m_Input.vbIsPass.size(); i++)
	{
		if(m_Input.vbIsPass[i])		fout << m_Input.vsRaw[i] << endl;
	}
	fout.close();

	ReportBamInfo();

	return true;
}

bool CINFO::ReportBamInfo()
{
	ofstream fout;
	fout.open((m_sOutputFile+".baminfo").c_str());
	
	fout << "{" << endl;
	fout << "\t\"sampleName\" : \"" << m_sBamFile.substr(m_sBamFile.find_last_of('/')+1) << "\","<< endl;
	fout << "\t\"bamInfo\" : [" << endl;

	for(unsigned int i=0; i<m_Input.vsChr.size(); i++)
	{
		if(!m_Input.vbIsPass[i])	continue;
		fout << "\t\t{" << endl;
		fout << "\t\t\t\"chr\" : \"" 		<< m_Input.vsChr[i] << "\"," << endl;
		fout << "\t\t\t\"pos\" : " 			<< m_Input.vnPos[i] << "," << endl;
		fout << "\t\t\t\"ref\" : \"" 		<< m_Input.vsRef[i] << "\"," << endl;
		fout << "\t\t\t\"alt\" : \"" 		<< m_Input.vsAlt[i] << "\"," << endl;
		fout << "\t\t\t\"vaf\" : "	 		<< m_Input.vfVaf[i] << "," << endl;
		fout << "\t\t\t\"strandBias\" : "	<< m_Input.vfStrandBias[i] << "," << endl;
		fout << "\t\t\t\"readLen\" : [";
		for(unsigned int j=0; j<m_Input.v2nReadLen[i].size()-1; j++)
			fout << m_Input.v2nReadLen[i][j] << ",";
		fout << m_Input.v2nReadLen[i][m_Input.v2nReadLen[i].size()-1] << "]," << endl;
		fout << "\t\t\t\"mapQ\" : [";
		for(unsigned int j=0; j<m_Input.v2nMapQ[i].size()-1; j++)
			fout << (int)m_Input.v2nMapQ[i][j] << ",";
		fout << (int)m_Input.v2nMapQ[i][m_Input.v2nMapQ[i].size()-1] << "]," << endl;
		fout << "\t\t\t\"baseQ\" : [";
		for(unsigned int j=0; j<m_Input.v2nBaseQ[i].size()-1; j++)
			fout << (int)m_Input.v2nBaseQ[i][j] << ",";
		fout << (int)m_Input.v2nBaseQ[i][m_Input.v2nBaseQ[i].size()-1] << "]" << endl;
		if(i==m_Input.vsChr.size()-1)
			fout << "\t\t}" << endl;
		else
			fout << "\t\t}," << endl;
	}
	fout << "\t]" << endl;
	fout << "}" << endl;

	fout.close();

	return true;
}

