#include "../header/varinfo.h"
#include "../header/bio_files.h"

using namespace std;
pthread_mutex_t mutex_temp = PTHREAD_MUTEX_INITIALIZER;

/*

- To Do:
  -- get read info for multiple loci
  
- Changelog:
  -- Date Mar-30-2021 SeokCholHong shulkhorn@gmail.com
     --- GetRead

*/

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
			string sSeq = "";
			for(int j=nPos; j<=nPos+1; j++){
				int nLoc = j - (b->core.pos+1);
				if(0<=nLoc){
					if(GetRead(b, nLoc)=="") break;
					sSeq+=GetRead(b, nLoc);
				}
			}
			cout << sSeq << endl;
		}


		// determine whether SNPs or Del-Ins


	}


	bam_close(finBam);

	return false;
}

string CINFO::GetRead(bam1_t *b, int nPos)
{
	string sSeq="";
	//cigar
	uint32_t nTemp=0;
	int l,nPl,nPreAlgn;
	for(l=b->core.l_qname, nPl=1, nPreAlgn = 0;
		l<(b->core.l_qname + (b->core.n_cigar*4));
		l++, nPl++){
		int nChk = nPl%4;
		if(nChk==0){
			nTemp = (b->data[l]<<24) | nTemp;

			int nAlgn = (nTemp>>4);
			int op = nTemp&0xf;
			if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF){
				nPreAlgn += nAlgn;
				if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CINS){
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CDEL){
				if(nPreAlgn<=nPos && nPos<nPreAlgn+nAlgn){
					return sSeq;
				}
				else if(nPos>=nPreAlgn+nAlgn){
					nPos -= nAlgn;
				}
				else if(nPos<nPreAlgn)	break;
			}
			else if(op == BAM_CSOFT_CLIP || op == BAM_CPAD){
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CREF_SKIP){
				if(nPos<nPreAlgn)	break;
				nPos += nAlgn;
				nPreAlgn += nAlgn;
			}
			else if(op == BAM_CHARD_CLIP){
			}
		}
		else if(nChk==1)	nTemp = b->data[l];
		else if(nChk==2)	nTemp = (b->data[l]<<8) | nTemp;
		else if(nChk==3)	nTemp = (b->data[l]<<16) | nTemp;
	}
	if(nPos>=b->core.l_qseq)	return sSeq;

	unsigned char ucBase = bam1_seqi(bam1_seq(b), nPos);
	if     (ucBase == 0x01)	sSeq="A";
	else if(ucBase == 0x02)	sSeq="C";
	else if(ucBase == 0x04)	sSeq="G";
	else if(ucBase == 0x08)	sSeq="T";
	else if(ucBase == 0x0F)	sSeq="N";
	else
		throw std::logic_error("ERROR: sequence error. It deosn't belong to ACGTN (in function GetRead())");

	return sSeq;
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




