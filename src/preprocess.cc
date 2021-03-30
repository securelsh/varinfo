#include "../header/varinfo.h"
#include "../header/bio_files.h"

using namespace std;
pthread_mutex_t mutex_temp = PTHREAD_MUTEX_INITIALIZER;


bool CINFO::ModVariant()
{
	for(int i=0; i<this->m_Input.vbIsVcf.size()-1; i++)
	{
		int nIdx = i;
		while(this->m_Input.vnPos[nIdx] + 1 == this->m_Input.vnPos[nIdx+1]
		&& this->m_Input.vsRef[nIdx].size() == 1 && this->m_Input.vsAlt[nIdx].size() == 1
		&& this->m_Input.vsRef[nIdx+1].size() == 1 && this->m_Input.vsAlt[nIdx+1].size() == 1)
		{
			nIdx++;
			if(nIdx == this->m_Input.vbIsVcf.size()-1)	break;
		}


		if(i != nIdx)			
		{
			this->DelIns(i, nIdx);
			i = nIdx;

			break;
		}


	}
	return false;
}


bool CINFO::DelIns(int nSIdx, int nEIdx)
{
	bamFile finBam;
	finBam = bam_open(m_sBamFile.c_str(), "r");
	if(finBam == NULL)	throw	std::logic_error("ERROR: Can not open .bam file");

	bam_header_t *bamHeader;
	bamHeader = bam_header_read(finBam);
	if(!bamHeader)		throw std::logic_error("ERROR: Fail to read the header of a bam file");

	bam_index_t *bamIndex;
	bamIndex = bam_index_load(m_sBamFile.c_str());
	if(!bamIndex)		throw std::logic_error("ERROR: Fail to read .bai file");

	string sChr = this->m_Input.vsChr[nSIdx];
	int nChr = ConvertChrToTid(sChr, bamHeader);
	int nPos = this->m_Input.vnPos[nSIdx];
	if(nChr == -1)
	{
		cout << sChr << "(" << nChr << "):" << nPos << endl;
		throw std::logic_error("Not available Chr");
	}



	for(int i=nSIdx; i<nEIdx; i++)
	{
		// extract reads
		bam_iter_t bamIter;
		bam1_t *b;
		b = bam_init1();
		bamIter = bam_iter_query(bamIndex, nChr, nPos-1, nPos);
		int nRet;
		while((nRet = bam_iter_read(finBam, bamIter, b)) >= 0)
		{
			GetRead(b, nPos, nPos+1);
			break;

		}


		// determine whether SNPs or Del-Ins


	}


	bam_close(finBam);

	return false;
}

bool CINFO::GetRead(bam1_t *b, int nSPos, nEPos)
{
	//cigar
	uint32_t nTemp = 0;
	int l, nPl, nPreAlign;
	for(l=b->core.l_qname, nPl=1, nPreAlign=0;
		l<(b->core.l_qname + (b->core.n_cigar*4));
		l++, nPl++)
	{
		int nChk = nPl%4;
		if(nChk==0){
			nTemp = (b->data[l]<<24) | nTemp;



		}
			
	}





	return false;
}



CINFO::CINFO(string sBamFile, string sInputFile, string sOutputFile, bool bIsMod, int nCntThread, double dVaf, bool bIsDebug)
{
	m_sBamFile = sBamFile;
	m_sInputFile = sInputFile;
	m_sOutputFile = sOutputFile;
	m_bIsMod = bIsMod;
	m_nCntThread = nCntThread;
	m_dVaf = dVaf;
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
			Parsing(vsWord, sLine, "\t");
			vector<string> vsDepth;
			Parsing(vsDepth, vsWord[1], " ");

			int nRefDp = strtol(vsDepth[0].c_str(), &pcEnd, 10);
			int nAltDp = strtol(vsDepth[1].c_str(), &pcEnd, 10);

			if((double)nAltDp/(double)(nRefDp+nAltDp) < this->m_dVaf)	continue;

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
			this->m_Input.vsRaw.push_back(sLine);
			this->m_Input.vbIsPass.push_back(true);
		}
	}
	
	this->ViewStatus(1,1,"Read Input File");
	cout << endl;


	return false;
}




