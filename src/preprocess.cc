#include "../header/varinfo.h"
#include "../header/bio_files.h"

using namespace std;
pthread_mutex_t mutex_temp = PTHREAD_MUTEX_INITIALIZER;


bool CINFO::ModVariant()
{

	return false;
}



CINFO::CINFO(string sBamFile, string sInputFile, string sOutputFile, bool bIsMod, int nCntThread, bool bIsDebug)
{
	m_sBamFile = sBamFile;
	m_sInputFile = sInputFile;
	m_sOutputFile = sOutputFile;
	m_bIsMod = bIsMod;
	m_nCntThread = nCntThread;
	m_bIsDebug = bIsDebug;
}


bool CINFO::ReadInput()
{
	//read vcf file
	this->SetStartTime();
	this->ViewStatus(0,1,"Read Input File");
	
	ifstream fin;
	fin.open(m_sInputFile.c_str());
	if(!fin.is_open())		throw std::logic_error("Input file does not exists");

	string sLine;
	char *pcEnd;
	vector<string> vsWord;
	while(getline(fin, sLine) && sLine.size() > 0)
	{
		if(sLine[0] == '#')			this->m_Input.vsHeader.push_back(sLine);
		else
		{
			this->m_Input.vsRaw.push_back(sLine);
			Parsing(vsWord, sLine, "\t");
			
			if(vsWord[0].substr(0,4) == "chr ")			//Adiscan format
			{
				this->m_Input.vbIsVcf.push_back(false);
				vector<string> vsSplit;
				Parsing(vsSplit, vsWord[0].substr(4), "_");
				this->m_Input.vsChr.push_back(vsSplit[0]);
				this->m_Input.vnPos.push_back(strtol(vsSplit[1].c_str(), &pcEnd, 10));
				
				Parsing(vsSplit, vsWord[2], " ");
				this->m_Input.vsRef.push_back(vsSplit[0]);
				this->m_Input.vsAlt.push_back(vsSplit[1]);
			}
			else										//Vcf format
			{
				this->m_Input.vbIsVcf.push_back(true);
				this->m_Input.vsChr.push_back(vsWord[0]);
				this->m_Input.vnPos.push_back(strtol(vsWord[1].c_str(), &pcEnd, 10));
				this->m_Input.vsRef.push_back(vsWord[3]);
				this->m_Input.vsAlt.push_back(vsWord[4]);
	
			}
		}
	}
	
	this->ViewStatus(1,1,"Read Input File");
	cout << endl;


	return false;
}




